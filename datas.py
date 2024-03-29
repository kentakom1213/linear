import linear
from fractions import Fraction
from functools import reduce
import math

### sqrt class
class Sqrt:
    def __init__(self, dec, factors=None):
        if dec == 0:
            return 0

        self.factors = factors or self.factoring(dec)
        self.rational, self.irrational = self.separate()
    
    def __float__(self):
        return self.rational * math.sqrt(self.irrational)

    def __mul__(self, other):
        if isinstance(other, int):
            other = Sqrt(other)
        merged_factors = self.merge_factors(self.factors, other.factors)
        return Sqrt(None, merged_factors)
    
    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, num):
        new_factors = [(n, i*num) for n, i in self.factors]
        return Sqrt(None, new_factors)
    
    def __repr__(self):
        if self.irrational == 1:
            return f"{self.rational}"
        elif self.rational == 1:
            return f"Sqrt({self.irrational})"
        else:
            return f"{self.rational}*Sqrt({self.irrational})"
    
    def separate(self):
        rational, irrational = 1, 1
        for n, cnt in self.factors:
            rational *= pow(n, cnt // 2)
            irrational *= pow(n, cnt % 2)
        return rational, irrational

    @staticmethod
    def factoring(n):
        result = []
        for i in range(2, int(n ** 0.5 + 1)):
            if n % i != 0: continue
            counter = 0
            while n % i == 0:
                n //= i
                counter += 1
            result.append((i, counter))
        if n != 1:
            result.append((n, 1))
        return result
    
    @staticmethod
    def merge_factors(f1, f2):
        dic = {n:cnt for n, cnt in f1}
        for n, cnt in f2:
            if n in dic:
                dic[n] += cnt
            else:
                dic[n] = cnt
        return list(dic.items())


### vector class
class Vector(list):
    def __init__(self, list1D, frac=True):
        if frac:
            self.vector = [Fraction(v) for v in list1D]
        else:
            self.vector = list1D
        self.dimention = len(list1D)
    
    def is_scalar(func):
        def wrapper(inst, val):
            if isinstance(val, (int, float, Fraction, Sqrt)):
                return func(inst, val)
            else:
                raise ValueError
        return wrapper
    
    def norm(self):
        """ノルムを求める"""
        return Sqrt(self.dot(self))
    
    def __abs__(self):
        """abs()関数"""
        return self.norm()
    
    def __add__(self, other):
        return Vector([x + y for x, y in zip(self.vector, other.vector)])

    def __sub__(self, other):
        return Vector([x - y for x, y in zip(self.vector, other.vector)])

    @is_scalar
    def __mul__(self, scalar):
        return Vector([x * scalar for x in self.vector])
    
    @is_scalar
    def __rmul__(self, scalar):
        return Vector([x * scalar for x in self.vector])
    
    @is_scalar
    def __truediv__(self, scalar):
        if scalar == 0:
            raise ZeroDivisionError()
        else:
            return Vector([x / scalar for x in self.vector])
    
    def dot(self, other):
        return sum(i*j for i, j in zip(self, other))
    
    def __iter__(self):
        yield from self.vector
    
    def __getitem__(self, index):
        if isinstance(index, int):
            return self.vector[index]
        else:
            raise ValueError
    
    def copy(self):
        return Vector(self.vector.copy())
    
    def __str__(self):
        list2str = lambda a, b: f"{a}, {b}"
        str_vec = "[" + reduce(list2str, self.vector) + "]"
        return str_vec

    def __repr__(self):
        list2str = lambda a, b: f"{a}, {b}"
        str_vec = "Vector([" + reduce(list2str, self.vector) + "])"
        return str_vec
    
    def to_list(self):
        return [int(v) for v in self.vector]


### matrix class
class Matrix():
    def __init__(self, arr, frac=True):
        if frac:
            self.matrix = [[Fraction(v) for v in row] for row in arr]
        else:
            self.matrix = arr
        self.shape = (len(arr), len(arr[0]))

        # check the shape of arr
        for i in range(self.shape[0]):
            if len(self.matrix[i]) != self.shape[1]:
                raise ValueError("Incorrect Input")
    
    def filter_matrix(func):
        def wrapper(inst, other):
            if isinstance(other, Matrix):
                return func(inst, other)
            else:
                raise ValueError("Value should be matrix.")
        return wrapper

    @filter_matrix
    def __add__(self, other):
        res = [[0] * self.shape[1] for _ in range(self.shape[0])]
        for i in range(self.shape[0]):
            for j in range(self.shape[1]):
                res[i][j] = self.matrix[i][j] + other.matrix[i][j]
        return Matrix(res)

    @filter_matrix
    def __sub__(self, other):
        res = [[0] * self.shape[1] for _ in range(self.shape[0])]
        for i in range(self.shape[0]):
            for j in range(self.shape[1]):
                res[i][j] = self.matrix[i][j] - other.matrix[i][j]
        return Matrix(res)
    
    @filter_matrix
    def __matmul__(self, other):
        """@ 演算子でMatrix同士の積を求める"""
        if self.shape[1] != other.shape[0]:
            raise ValueError("The matrix shape is inappropriate.")
        else:
            res = [[0] * other.shape[1] for _ in range(self.shape[0])]
            for i in range(self.shape[0]):
                for j in range(other.shape[1]):
                    res[i][j] = self[i, :].dot(other[:, j])
            return Matrix(res)
    
    def __pow__(self, num):
        if isinstance(num, int):
            res = self[:, :]
            for _ in range(num-1):
                res @= res
            return res
        else:
            ValueError
    
    def __mul__(self, scalar):
        res = [[0] * self.shape[1] for _ in range(self.shape[0])]
        for i in range(self.shape[0]):
            for j in range(self.shape[1]):
                res[i][j] = self.matrix[i][j] * scalar
        return Matrix(res)
    
    def __rmul__(self, scalar):
        return self.__mul__(scalar)
    
    def __truediv__(self, scalar):
        res = [[0] * self.shape[1] for _ in range(self.shape[0])]
        for i in range(self.shape[0]):
            for j in range(self.shape[1]):
                res[i][j] = self.matrix[i][j] / scalar
        return Matrix(res)
    
    def transpose(self):
        """転置行列を求める"""
        res = [[0] * self.shape[0] for _ in range(self.shape[1])]
        for i in range(self.shape[0]):
            for j in range(self.shape[1]):
                res[j][i] = self.matrix[i][j]
        return Matrix(res)
    
    def row_reduction(self):
        res = linear.row_reduction(self.matrix)
        return Matrix([row.to_list() for row in res])
    
    def inverse(self):
        res = [row.to_list() for row in linear.find_inverse_matrix(self.matrix)]
        return Matrix(res)
    
    def cofactor(self, i, j):
        """余因子行列を求める"""
        arr = self.matrix
        res = [[el for el in r[:j]+r[j+1:]] for r in arr[:i]+arr[i+1:]]
        return Matrix(res)
    
    def det(self):
        """行列式を余因子展開によって求める"""
        arr = self.matrix
        N = len(arr)
        if N == 1:
            return arr[0][0]
        else:
            det_arr = 0
            for i in range(N):
                det_arr += (-1)**i * arr[0][i] * self.cofactor(0, i).det()
            return det_arr
    
    def __iter__(self):
        yield from map(Vector, self.matrix)
    
    def __getitem__(self, item):
        rslice, cslice = item
        if isinstance(rslice, int) and isinstance(cslice, int):
            return self.matrix[rslice][cslice]
        elif isinstance(rslice, int):
            return Vector(self.matrix[rslice])
        elif isinstance(cslice, int):
            res = [row[cslice] for row in self.matrix]
            return Vector(res)

        rstart, rstop = (rslice.start, rslice.stop)
        cstart, cstop = (cslice.start, cslice.stop)
        if rstart == rstop != None or cstart == cstop != None:
            raise ValueError("The interval is invalid")

        new_arr = [row[cstart:cstop] for row in self.matrix[rstart:rstop]]
        return Matrix(new_arr)

    def __repr__(self):
        res = "Matrix(["
        for i, row in enumerate(self.matrix):
            if i > 0:
                res += "        "
            res += str(Vector(row))
            if i < self.shape[0]-1:
                res += ",\n"
            else:
                res += "])"
        return res
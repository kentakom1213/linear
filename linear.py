from fractions import Fraction
from datas import *

### row_reduction ###
def row_reduction(array, show=True):
    array = [Vector([Fraction(el) for el in list]) for list in array]
    dimention = len(array)
    
    for i, _ in enumerate(array):
        if show:
            print(*map(str, array), sep=", ")
        if array[i][i] == 0:
            array[i], array[-1] = array[-1], array[i]
        #print(f"i = {i}") ##### debug #####
        try:
            array[i] = _make_one(array[i], i)
        except ZeroDivisionError:
            return
        rest = list(range(i)) + list((range(i+1, dimention)))
        for j in rest:
            array[j] = _make_zero(array[j], array[i], i)
    
    result = array #[list(map(float, vec.to_list())) for vec in array]
    return result


def find_inverse_matrix(array):
    degree = len(array)
    identity = get_identity(degree)
    extended_array = [row1 + row2 for row1, row2 in zip(array, identity)]
    reducted = row_reduction(extended_array, show=False)
    if reducted == None:
        print("This matrix is not regular.")
        return
    inverse = [Vector(row.to_list()[degree:]) for row in reducted]
    return inverse

def _make_one(row: "list", row_number: "int") -> "list":
    base = row[row_number]
    return row / base

def _make_zero(target_line, criterion_line, row_number):
    """
    In : [1, 3, 2], [0, 1, 5], 1
    -> [1, 3, 2] += [0, 1, 5] * (-3)
    Out: [1, 0, -13]"""
    multiply_number = - ( target_line[row_number] / criterion_line[row_number] )
    return target_line + criterion_line * multiply_number

def solve(array):
    neg_last_item = lambda x: x[:-1] + [-x[-1]]
    # array = list( map(neg_last_item, array) )
    solved = row_reduction(array, show=False)
    return Vector([x[-1] for x in solved])


### array ###
def get_zeros(degree):
    return [[0]* degree for _ in range(degree)]

def get_identity(degree):
    identity = get_zeros(degree)
    for i in range(degree):
        identity[i][i] = 1
    return identity


### det
cofactor = lambda a, i, j: [[el for el in r[:j]+r[j+1:]] for r in a[:i]+a[i+1:]]

def det(arr):
    N = len(arr)
    if N == 1:
        return arr[0][0]
    else:
        det_arr = 0
        for i in range(N):
            det_arr += (-1)**i * arr[0][i] * det(cofactor(arr, 0, i))
        return det_arr


### dot
def dot(a, b):
    return sum(i*j for i, j in zip(a, b))


if __name__ == "__main__":
    A = [[1, 2, 3, 4, 5, -1],
         [2, 2, 2, 2, 4, -3],
         [0, 4, 1, 1, 3, -1],
         [2, -1, 2, 1, 2, -4],
         [3, -2, 2, 1, 1, -2]]

    print(solve(A))

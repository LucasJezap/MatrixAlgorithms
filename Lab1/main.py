import random
import numpy as np
from datetime import datetime

def create_matrix(M, a, b, generate):
    """
    Function which reads matrix values
    :param M: matrix
    :param a: number of rows
    :param b: number of columns
    :param generate: boolean value
    :return: matrix with values
    """
    if not generate == 'y':
        print("Please enter consecutive matrix rows "
          "(the size of the matrix is " + str(a) + "x" + str(b) + ")")
    for i in range(a):
        if generate == 'y':
            M[i] = [random.randint(0,100) - random.randint(0,100) for i in range(b)]
        else:
            M[i] = list(map(lambda x : int(x), input().split()))

    if generate == 'y':
        print("Generated matrix")
        print(np.matrix(M))
    return M

def dot_product(row, column, k):
    """
    function returning the dot product of row and column
    :param row: given row
    :param column: given column
    :return: the dot product of row and column
    """
    res = 0
    for i in range(k):
        res += row[0,i] * column[i,0]
    return res

def matrix_multiplication_jip(A, B, m, k, n):
    """
    Dense matrix multiplication of matrixes A and B (version 3a in the script)
    :param A: matrix A
    :param B: matrix B
    :param m: number of rows of A
    :param k: number of columns of A and number of rows of B
    :param n: number of columns of B
    :return: result matrix C
    """
    C = [[0] * n for i in range(m)]
    start = datetime.now()
    for j in range(n):
        for i in range(m):
            C[i][j] += dot_product(A[i], B[:,j], k)

    end = datetime.now()
    print("Time3 = " + str(end - start))
    return C


if __name__ == '__main__':
    generate = input("Do you want to generate the matrixes? (y/n) \n")
    if generate != 'y' and generate != 'n':
        print("Wrong answer")
        exit()
    print("Matrix A:")
    m = int(input("Please enter number of rows: "))
    k = int(input("Please enter number of columns: "))
    A = np.matrix([ [0]*k for  i in range(m)])
    create_matrix(A, m, k, generate)
    print("Matrix B:")
    tmp = int(input("Please enter number of rows: "))
    if tmp != k:
        print("Number of columns in A and number of rows in B must be the same!")
        exit()
    n = int(input("Please enter number of rows: "))
    B = np.matrix([[0] * n for i in range(k)])
    create_matrix(B, k, n, generate)
    start = datetime.now()
    C = matrix_multiplication_jip(A, B, m, k, n)
    end = datetime.now()
    print("\n\n\nMatrix C:")
    print(np.matrix(C))
    print("Numpy calculated C:")
    start2 = datetime.now()
    np.matmul(np.array(A), np.array(B))
    end2 = datetime.now()
    print("Time = " + str(end-start))
    print("Time2 = " + str(end2 - start2))
    print(str((end - start) / (end2-start2)))


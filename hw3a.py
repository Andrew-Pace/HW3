#region imports
import Gauss_Seidel as GS
import DoolittleMethod as dm
import random
import math
#endregion

#region functions
def Cholesky(Aaug):
    """
    This function finds the solution to a matrix equation Ax=b by the Colesky method
    :param Aaug: An augmented matrix
    :return: the solution vector x, L and Ltrans as a tuple
    """
    #step 1:  split the Aaug into A and b (see separateAugmented in Gauss_Seidel.py)
    A, b = GS.separateAugmented(Aaug) # using gauss_seidel.py function to seperate augmented matrix
    n = len(A)
    for i in range(0, len(Aaug) - 1):
        A.append(Aaug[i])

    #step 2:  factor into L and Ltrans according to Cholesky formulae
    N = len(A) - 1 # getting length - 1
    lower = [[0 for x in range(N + 1)]
             for y in range(N + 1)]

    # Decomposing A matrix into Lower Triangular
    for i in range(n):
        for j in range(i + 1):
            sum1 = 0
            if (j == i):
                for k in range(j):
                    sum1 += pow(lower[j][k], 2)
                lower[j][j] = int(math.sqrt(A[j][j] - sum1))
            else:
                for k in range(j):
                    sum1 += (lower[i][k] * lower[j][k])
                if (lower[j][j] > 0):
                    lower[i][j] = int((A[i][j] - sum1) / lower[j][j])

    #step 3:  use backsolving to find x (see methods in DoolittleMethod.py)
    # reusing code from DoolittleMethod.Backsolve to backsolve Ly=b
    ynRows=len(b) # marking variables as y for first step
    ynCols=ynRows
    y=[0]*ynRows
    UT=False # using lower tri
    if UT: # didnt need to leave in for current matrices but may increase versatility in the future?
        for nR in range(ynRows-1,-1,-1):
            s=0
            for nC in range(nR+1,ynRows):
                s+=lower[nR][nC]*y[nC]
            y[nR]=1/lower[nR][nR]*(b[nR]-s)
    else:
        for nR in range(ynRows):
            s=0
            for nC in range(nR):
                s+=lower[nR][nC]*y[nC]
            y[nR]=1/lower[nR][nR]*(b[nR]-s)

    # Backsolving for Ltp*x=y, code from doolittles method.backsolve
    nRows=len(b)  # variables for x for second step
    nCols=nRows
    x=[0]*nRows
    UT=True # setting upper tri to true for transposed calculations
    Ltp = Transpose(lower) # using transpose function for lower
    if UT:
        for nR in range(nRows-1,-1,-1):
            s=0
            for nC in range(nR+1,nRows):
                s+=Ltp[nR][nC]*x[nC]
            x[nR]=1/Ltp[nR][nR]*(y[nR]-s)
    else: # didnt need to leave in for current matrices but may increase versatility in the future?
        for nR in range(nRows):
            s=0
            for nC in range(nR):
                s+=Ltp[nR][nC]*x[nC]
            x[nR]=1/Ltp[nR][nR]*(y[nR]-s)

    #step 4:  return (x,L,Ltrans)
    return x, lower, Ltp

def SymPosDef(A):
    """
    This function first finds the transpose of A and then compares all elements of A to Atrans.
    If I pass that test, I create a vector x with random numbers and perform xtrans*A*x to see if>0.
    :param A: a nxn matrix
    :return: True if symmetric, positive definite
    """
    #step 1:  recall that a transpose has elements such that Atrans[i][j] = A[j][i] see page 267 in MAE3013 text
    tp = Transpose(A)
    if tp == A:     #step 2:  check that all elements of A and Atrans are the same. if fail->return false
        x = []
        n = len(A[0])
        for i in range(n):
            x.append(random.uniform(-1, 1))     #step 3:  produce a vector x of length n filled with random floats between -1 and +1
        result = []
        for m in range(0, len(A)):    #step 4:  compute xtrans*A*x
            rows = []
            for i in range(0, 1):
                columns = 0
                for j in range(0, len(A)):
                    columns += A[m][j] * x[j]
                rows.append(columns)
            result.append(rows)
        sum = 0
        for m in range(0, n):
            sum += result[m][0] * x[m]
        if sum > 0:  #step 5:  if step 4 > 0 return true else return false
            return True
        else:
            return False

    else:
        return False
        print("x:", x)
        print("Ax:", result)
        print("sum:", sum)

def Transpose(A):
    """
    This function finds the transpose of a square matrix
    :param A: an nxn matrix
    :return: the transpose of A
    """
    N = len(A[0]) # defining size of matrix
    M = len(A) # defining size of matrix
    B = [[0 for x in range(M)] for y in range(N)] # setting variable
    for i in range(N):
        for j in range(M):
            B[i][j] = A[j][i] # transposing

    # print("Result matrix is") #used to print transpose for testing purposes
    # for i in range(N):
    #     for j in range(M):
    #         print(B[i][j], " ", end='')
    #     print()
    return B

def main():
    """
    Step 1:  I need to first define the matrices given in part a) of HW3_2024.
    Step 2:  pass a matrix to SymPosDef to tell if it is symmetric, positive definite
    Step 3:  based on result of Step 2, use either the Doolittle or Cholesky method to solve
    Steo 4:  check my answer by multiplying A*x to see if I get b
    Step 4:  print the solution vector and which method was used to the cli
    """
    Aaug = [[1,-1,3,2,15] , [-1,5,-5,-2,-35], [3,-5,19,3,94], [2,-2,3,21,1]] # defining the augmented matrix
    A, b = GS.separateAugmented(Aaug) # seperating Aaug
    if SymPosDef(A) == True: # checking if it passes SymPosDef

        # Making it look pretty
        print("Cholesky Method Used")
        for i in range(len(A[0])):
            for j in range(len(A)):
                print(A[i][j], " ",end='')
                if j == len(A)-1: print(b[i], end='')
            print()
        print("Solution Vector: ",Cholesky(Aaug)[0])
        print()

    else: # using doolittle if does not pass SymPosDef

        # Making it look pretty
        print("Doolittle Method Used")
        for i in range(len(A[0])):
            for j in range(len(A)):
                print(A[i][j], " ",end='')
                if j == len(A)-1: print(b[i], end='')
            print()
        print("Solution Vector: ",dm.Doolittle(Aaug))
        print()


    Aaug = [[4,2,4,0,20] , [2,2,3,2,36], [4,3,6,3,60], [0,2,3,9,122]]  # defining the augmented matrix
    A, b = GS.separateAugmented(Aaug) # seperating Aaug
    if SymPosDef(A) == True: # checking if it passes SymPosDef

        # Making it look pretty
        print("Cholesky Method Used")
        for i in range(len(A[0])):
            for j in range(len(A)):
                print(A[i][j], " ",end='')
                if j == len(A)-1: print(b[i], end='')
            print()
        print("Solution Vector: ",Cholesky(Aaug)[0])
        print()
    else:

        # Making it look pretty
        print("Doolittle Method Used")
        for i in range(len(A[0])):
            for j in range(len(A)):
                print(A[i][j], " ",end='')
                if j == len(A)-1: print(b[i], end='')
            print()
        print("Solution Vector: ",dm.Doolittle(Aaug))
        print()

    # for testing doolittles method
    # Aaug = [[1,1,-1,4] , [1,-2,3,-6], [2,3,1,7]]  # defining the augmented matrix
    # A, b = GS.separateAugmented(Aaug)
    # if SymPosDef(A) == True: # checking if it passes SymPosDef
    #     # Making it look pretty
    #     print("Cholesky Method Used")
    #     for i in range(len(A[0])):
    #         for j in range(len(A)):
    #             print(A[i][j], " ",end='')
    #             if j == len(A)-1: print(b[i], end='')
    #         print()
    #     print("Solution Vector: ",Cholesky(Aaug)[0])
    # else:
    #     # Making it look pretty
    #     print("Doolittle Method Used")
    #     for i in range(len(A[0])):
    #         for j in range(len(A)):
    #             print(A[i][j], " ",end='')
    #             if j == len(A)-1: print(b[i], end='')
    #         print()
    #     print("Solution Vector: ",dm.Doolittle(Aaug))

if __name__ == "__main__":
    main()
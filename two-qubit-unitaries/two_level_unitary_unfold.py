import numpy as np

'''

For a NxN matrix, the following function takes the matrix as input and makes the
(0,0) element 1 and all other elements zero in the zeroth row and the zeroth column.
The function returns a pair where the first elementi s list of the two level unitary
gates U such that U_(N-1)*...*U_1 = (matrix with 1 in (0,0) position) and the second
element is the product U_(N-1)*...*U_1. It follows the algorithm outlined in
Neilsen and Chuang (page 189).

'''
def ReductionN(Matrix):

    Two_Level_Universal_Gate = []
    
    # Extract the size of the input matrix.
    Matrix_dim = int(np.sqrt(Matrix.size))

    for i in range(1,Matrix_dim):

        if Matrix[i,0] == 0:

            U = np.identity(Matrix_dim, dtype = complex)

            if i == 1:

                pass

            else:

                U.itemset((0,0), np.conjugate(Matrix[0,0]))

        else:

            a, b = Matrix[0,0], Matrix[i,0]

            u00 = np.conjugate(a) / np.sqrt(np.absolute(a)**2 + np.absolute(b)**2)
            u01 = np.conjugate(b) / np.sqrt(np.absolute(a)**2 + np.absolute(b)**2)
            u10 = b / np.sqrt(np.absolute(a)**2 + np.absolute(b)**2)
            u11 = - a / np.sqrt(np.absolute(a)**2 + np.absolute(b)**2)
    
            U = np.identity(Matrix_dim, dtype = complex)
            U.itemset((0,0), u00)
            U.itemset((0,i), u01)
            U.itemset((i,0), u10)
            U.itemset((i,i), u11)
        
        Two_Level_Universal_Gate.append(U)

        # Calculating the product of all the U matrices with G.

        Matrix = np.matmul(U, Matrix)

    return Two_Level_Universal_Gate, Matrix


'''

The following function takes a 3x3 matrix U as input and returns a list of matrices
U_1, U_2 and U_3 such that U_3*U_2*U_1*U = identity (3x3). Follows the algorithm outlined
in Neilsen and Chuang (page 189).

'''    

def Reduction3(Matrix):

    Two_Level_Universal_Gate = []
    
    # Size of the matrix is 3.
    Matrix_dim = 3
 
    for i in range(1,Matrix_dim):

        '''
            If the (1,0) element is 0, then U_1 = identity matrix.
        
        '''
        
        if Matrix[i,0] == 0:

            U = np.identity(Matrix_dim, dtype = complex)

            
            if i == 1:
                
                pass # For the (1,0) element U_1 is identity.
            
            else:
                
                U.itemset((0,0), np.conjugate(Matrix[0,0]))

        else:

            a, b = Matrix[0,0], Matrix[i,0]

            u00 = np.conjugate(a) / np.sqrt(np.absolute(a)**2 \
                + np.absolute(b)**2)
            u01 = np.conjugate(b) / np.sqrt(np.absolute(a)**2 \
                + np.absolute(b)**2)
            u10 = b / np.sqrt(np.absolute(a)**2 + np.absolute(b)**2)
            u11 = - a / np.sqrt(np.absolute(a)**2 + np.absolute(b)**2)
    
            U = np.identity(Matrix_dim, dtype = complex)
            U.itemset((0,0), u00)
            U.itemset((0,i), u01)
            U.itemset((i,0), u10)
            U.itemset((i,i), u11)
        
        Two_Level_Universal_Gate.append(U)

        # Calculating product of all the U matrices with G.
        Matrix = np.matmul(U, Matrix)

    # U3 is calculated sepatarely.
    Two_Level_Universal_Gate.append(np.conjugate(Matrix).transpose())

    return Two_Level_Universal_Gate

'''

The following function takes a matrix as input and returns a matrix with the
first row and the first column deleted.


'''

def Extract_subMatrix(Matrix):

    Matrix_size = int(np.sqrt(Matrix.size))
    subMatrix = Matrix[1:Matrix_size, 1:Matrix_size]

    return  subMatrix


'''

The following function takes a matrix of size NxN and an integer K>N and returns a matrix
whose size is KxK, with extra row and columns same as that of identity matrix of size KxK.

'''
def Expand_Matrix(Matrix, Required_Size):

    # Intital size of the input matrix.
    Initial_Size = int(np.sqrt(Matrix.size))

    # Creating an initial identity matrix.
    M = np.identity(Required_Size, dtype = complex)

    # Size difference between the two matrics.
    Sz = Required_Size - Initial_Size

    # The elements from the old matrix is added to M, other row and column are part of
    # identity matrix, this gurantees that the matrix is a two level unitary matrix.

    for i in range(Sz, Required_Size):
        
        for j in range(Sz, Required_Size):
            
            M.itemset((i,j), Matrix[i-Sz,j-Sz])

    return M

def Reduce(Matrix):
    N = np.int(np.sqrt(Matrix.size))
    # This list will contain all the two level unitary gates.
    Unitary_gates = []

    Last_Reduced_Matrix = []

    for i in range(N, 2, -1): # Loop from N to 3.
    
    
        # 3x3 matrices are treated separetely with CReduction3(Matrix).
        if i == 3:
        
            M = Last_Reduced_Matrix[-1]

            for g in Reduction3(M):
            
                Unitary_gates.append(g)

        # First iteration, start with original matrix G.
        elif i == N:

            M = Matrix
            
            for g in ReductionN(M)[0]: # List of two level gates.
            
                Unitary_gates.append(g)

        
            # The reduced submatrix is put into the list to use in the next loop.
            Last_Reduced_Matrix.append(Extract_subMatrix(ReductionN(M)[1]))
            #print(Last_Reduced_Matrix)

        else:
        
            # Any other dimension is reduced with CReductionN(Matrix).
            M = Last_Reduced_Matrix[-1]

            for g in ReductionN(M)[0]:
            
                Unitary_gates.append(g)   

            # The reduced submatrix is put into the list to use in the next loop.
            Last_Reduced_Matrix.append(Extract_subMatrix(ReductionN(M)[1]))        
            #print(Last_Reduced_Matrix)
    return Unitary_gates


def Unfold(Matrix):

    N = np.int(np.sqrt(Matrix.size))

    Gates = []
    i = 1
    for g in Reduce(Matrix):

        # If the matrix is of size NxN, print the matrix.
        if int(np.sqrt(g.size)) == N :
            
            Gates.append(["U"+str(i), g])
            #print('U'+str(i)+'\n', g,'\n\n')

        else:

            # If the matrix is of lower size other than N, expand the matrix into a NxN
            # matrix.
            
            Gates.append(["U"+str(i), Expand_Matrix(g, N) ])
            #print('U'+str(i)+'\n', Expand_Matrix(g, N),'\n\n')

        i += 1
        
    return Gates                   
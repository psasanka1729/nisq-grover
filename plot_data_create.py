'''
This code creates a file : index, Energy of the bulk, Average entropy, std of entropy.
It takes many data files with same values of epsilon but different noise realizations
(different random numbers each time), and calculates the average entropy corresponding
to each eigenvalues and then put them into the file.

'''


'''
Input : One single file with same values of epsilon but different noise realizations.
Output : Text file with index, energy, entropy, std of entropy.

'''
import numpy as np

'''

                J = Index of Eigenvalue 
                P = Sorted Eigenvalue
                S = Average Entanglement Entropy
                

'''
N = 8

# If the bulk and the two special states comes closer than this, then they will be
# declared merged.
threshold = 1.e-6

'''
    J = index of eigenvalue.
    E = eigenvalue.
    S = entropy.

'''
J,E,S = np.loadtxt('sorted_entropy_0.000005.txt', delimiter = ',', unpack=True)

M = int(len(S)/2**N)

'''

 We will reshape the array S into a (M, 2**N) matrix. The i-th column of this matrix are
 the list of entropies cooresponding to i-th eigenvalue for different noise realizations.
 
 
 
 

'''
Sm = np.matrix(np.reshape(S,(M, 2**N)))

'''

We will reshape the array P into a (M, 2**N) matrix.

'''

Em = np.matrix(np.reshape(E,(M, 2**N)))

E_average = np.sum(Em, axis = 0)/M

def Array2List(Arr):
    Arr_l = Arr.tolist()
    l = []
    for i in Arr_l:
        l.append(i[0])
    return l

'''

The following function returns True if bulk and special states have merged,
otherwise returns False.

'''

def Bulk_Special_Merge(Matrix): # Matrix is Phi i.e. the list of SORTED eigenstates.
    
    for i in range(M): # for each phi
        
        x1 = Matrix[i:i+1].tolist()[0][0] # first element
        x2 = Matrix[i:i+1].tolist()[0][1] # second element
        
        x_2N = Matrix[i:i+1].tolist()[0][-1] # last element
        x_2N_1 = Matrix[i:i+1].tolist()[0][-2] # second last element
        
        if abs(x1-x2) < threshold or abs(x_2N - x_2N_1) < threshold: 
            
            return True # special states and bulk have merged
        
        else:
            
            return False # special states and bulk have not merged



if Bulk_Special_Merge(Em) == False: # Have not merged
        
    # We need to ignore the first and the last column (the two special states).
    
    Sm_bulk = Sm[:,1:2**N-1] # Sm_bulk has the bulk states.
    E_bulk = E_average[:,1:2**N-1]

    ind = [i for i in range(1,2**N-1)] # there are 2**N-2 eigenvalues.
    print('Have not Merged')
        
else:
    
    Sm_bulk = Sm
    E_bulk = E_average
    ind = [i for i in range(2**N)] # there are 2**N eigenvalues.
    print("Merged")

#sum_s = np.sum(Sm_bulk, axis = 0)
sum_s = np.sum(Sm_bulk, axis = 0) # sums the columns and returns an array of length 2**N.

'''

To calculate the average entropy for an eigenvalue, we divide the sum of all
 entropies for different noise 
(the array above we just calculated)by total number of noise realizations (M).

'''

average_s = sum_s/M # averaging over noise realization.

Sm_std = np.std(Sm_bulk, axis = 0)/np.sqrt(len(Sm_bulk)-1) # calculates the standard deviation of entropy array.

pf = open('PlotData.txt','w')
pf = open('PlotData.txt','a')

for i in range(len(ind)):

    pf.write(str(ind[i])+'\t'+ str(Array2List(E_bulk[:,i])[0]) +'\t' + str(Array2List(average_s[:,i])[0])+ 
        '\t' + str(Array2List(Sm_std[:,i])[0])+'\n') 


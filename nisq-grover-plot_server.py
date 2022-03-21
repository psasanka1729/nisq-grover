#!/usr/bin/env python
# coding: utf-8

# In[32]:

import sys
from qiskit import*
from qiskit import Aer
import qiskit.quantum_info as qi
import numpy as np
import re
#import csv
import time
#start = time.time()


# ### Target state $|w>$

# In[33]:

N_eps = 200      #int(sys.argv[1])
Target_state = '0000'      #sys.argv[2]



# ## Setting up the The Grover operator
# $$ G = U_w U_s$$
# Where 
# $$U_w = 2|w><w| - I$$ 
# $U_w$ is a matrix with the position of the target state 1 and all its diagnoal elements are -1 and rest all zero.

# In[34]:


# First we note the length of N.
N = len(Target_state)


## The operator U_s.
A = np.ones((2**N, 2**N))
U_s = (2/(2**N))*A - np.identity(2**N, dtype = complex)


## The operator U_w. This is neeed for the sign adjustment of Grover_reconstructed operator.
U_w = - np.identity(2 ** N, dtype=complex) 
Target_index = int(Target_state, 2)
U_w.itemset((Target_index, Target_index),1)


## G is the Grover operator.
G = np.matmul(U_w, U_s)




# The following list has all the gates in the format [name of the gate, angle, qubit].
l = []

file1 = open('gates_list.txt', 'r')
Lines = file1.readlines()
 

for line in Lines:
    l.append(line.strip())

gates_list = []

Rz_Number = 0
for i in range(len(l)):
    
    l_temp = []
    gate_name = l[i].split(',')[0]
    if gate_name == 'rz':
        Rz_Number +=1 
    gate_angle = l[i].split(',')[1]
    gate_qubit = l[i].split(',')[2]
    
    l_temp.append(gate_name)
    l_temp.append(gate_angle)
    l_temp.append(gate_qubit)

    gates_list.append(l_temp)
    


# ## The basis gates
# The following returns the matrix of the Hadamard, CNOT and RZ gate.

# ### Hadamard gate

# In[39]:


## The dimension of the matrix is fixed by the number of qubits.
def Hadamard(Qubit):
    
    ## Changing the simulator 
    backend = Aer.get_backend('unitary_simulator')

    ## The circuit without measurement
    circ = QuantumCircuit(N)
    circ.h(Qubit)

    ## job execution and getting the result as an object
    job = execute(circ, backend)
    result = job.result()

    ## get the unitary matrix from the result object
    return result.get_unitary(circ)


# ### CNOT gate

# In[40]:


## The dimension of the matrix is fixed by the number of qubits.
def CNOT(t,c):
    ## Changing the simulator 
    backend = Aer.get_backend('unitary_simulator')

    ## The circuit without measurement
    circ = QuantumCircuit(N)
    circ.cx(t,c)

    ## job execution and getting the result as an object
    job = execute(circ, backend)
    result = job.result()

    ## get the unitary matrix from the result object
    return result.get_unitary(circ) 


# ### RZ gate

# In[41]:


def Rz(Angle, Qubit):
    ## Changing the simulator 
    backend = Aer.get_backend('unitary_simulator')

    ## The circuit without measurement
    circ = QuantumCircuit(N)
    circ.rz(Angle, Qubit)

    ## job execution and getting the result as an object
    job = execute(circ, backend)
    result = job.result()

    ## get the unitary matrix from the result object
    return result.get_unitary(circ)     
# ## Adding noise to the Oracle

# ### Noise creation
# This following line creates an array of random numbers which length is equal to the number of Rz gates.
# $$\delta = [\delta_1, \delta_2, ..., \delta_{\text{Number of Rz}}]$$
# Where $\delta_i$ are random numbers between $-1$ and $1$.  



Rz_Noise = 2*(np.random.rand(Rz_Number)-0.5)


# In the following $\epsilon$ acts as a strength of noise. The array Rz_Noise has random numbers from -1 to 1. The following function multiplies $\epsilon$ to each random numbers in Rz_Noise and add that to each Rz gate of the oracle $U_w$. We will first reconstruct the Grover operator by first constructing the oracle $U_w$ from the transpile and then by multiplying $U_w$ and $U_s$ to get $G$.

# In[45]:


def Grover_reconstructed(epsilon):
    
    ## Initializing the oracle U_w as an identity matrix.
    
    Or = np.identity(2**N, dtype = complex) 

    ## In the following loop we multiply all the 1 and 2 qubit gates with (or without) noise.
    
    j = 0 # Index for the random noise list.
    
    for i in range(len(gates_list)): # l is the list with all the gates.
    
        if gates_list[i][0] == 'rz':
            
            Noise = np.random.rand(1)[0]
            Or = np.matmul(Or, Rz(float(gates_list[i][1])  + epsilon * Rz_Noise[j], int(gates_list[i][2])))
            
            j = j + 1
        
        elif gates_list[i][0] == 'h':
        
            Or = np.matmul(Or, Hadamard(int(gates_list[i][2])))
        
        else:
        
            Or = np.matmul(Or, CNOT(int(gates_list[i][1]), int(gates_list[i][2])))
     
    ## In the following we will fix the phase of the reconstructed Oracle. First we will make all the elements
    # 1 or -1.
    Or = Or/Or[0,0]
    
    ## The sign of the reconstructed Oracle should be same as that of original U_w.
    if np.sign(Or[0,0]) == np.sign(U_w[0,0]):
        
        pass # If the sign is same, then pass.
    
    else:
        
        Or = -Or # Otherwise change the sign.
    Gr = np.matmul(Or, U_s) ## The Grover operator G = U_w * U_s.
    
    return Gr



import numpy
import numpy.linalg

sigma_x=numpy.zeros((2,2),dtype=complex)
sigma_y=numpy.zeros((2,2),dtype=complex)
sigma_z=numpy.zeros((2,2),dtype=complex)
sigma_0=numpy.identity(2,dtype=complex)
sigma_x[0,1]=1.
sigma_x[1,0]=1.
sigma_y[0,1]=-1.j
sigma_y[1,0]=1.j
sigma_z[0,0]=1.
sigma_z[1,1]=-1.
sigma_plus=(sigma_x+1.j*sigma_y)/2.
sigma_minus=(sigma_x-1.j*sigma_y)/2.

def adjoint(psi):
    return psi.conjugate().transpose()

def psi_to_rho(psi):
    return numpy.outer(psi,psi.conjugate())

def exp_val(psi, op):
    return numpy.real(numpy.dot(adjoint(psi),op.dot(psi)))

def norm_sq(psi):
    return numpy.real(numpy.dot(adjoint(psi),psi))

def normalize(psi,tol=1e-9):
    ns=norm_sq(psi)**0.5
    if ns < tol:
        raise ValueError
    return psi/ns

def comm(a,b):
    return a.dot(b)-b.dot(a)

def anti_comm(a,b):
    return a.dot(b)+b.dot(a)

def is_herm(M,tol=1e-9):
    if M.shape[0]!=M.shape[1]:
        return False
    diff=M-adjoint(M)
    return max(numpy.abs(diff.flatten())) < tol

def is_unitary(M,tol=1e-9):
    if M.shape[0]!=M.shape[1]:
        return False
    diff=M.dot(adjoint(M))-numpy.identity((M.shape[0]))
    return max(numpy.abs(diff.flatten())) < tol

def eigu(U,tol=1e-9):
    (E_1,V_1)=numpy.linalg.eigh(U+adjoint(U))
    U_1=adjoint(V_1).dot(U).dot(V_1)
    H_1=adjoint(V_1).dot(U+adjoint(U)).dot(V_1)
    non_diag_lst=[]
    j=0
    while j < U_1.shape[0]:
        k=0
        while k < U_1.shape[0]:
            if j!=k and abs(U_1[j,k]) > tol:
                if j not in non_diag_lst:
                    non_diag_lst.append(j)
                if k not in non_diag_lst:
                    non_diag_lst.append(k)
            k+=1
        j+=1
    if len(non_diag_lst) > 0:
        non_diag_lst=numpy.sort(numpy.array(non_diag_lst))
        U_1_cut=U_1[non_diag_lst,:][:,non_diag_lst]
        (E_2_cut,V_2_cut)=numpy.linalg.eigh(1.j*(U_1_cut-adjoint(U_1_cut)))
        V_2=numpy.identity((U.shape[0]),dtype=V_2_cut.dtype)
        for j in range(len(non_diag_lst)):
            V_2[non_diag_lst[j],non_diag_lst]=V_2_cut[j,:]
        V_1=V_1.dot(V_2)
        U_1=adjoint(V_2).dot(U_1).dot(V_2)

    # Sort by phase
    U_1=numpy.diag(U_1)
    inds=numpy.argsort(numpy.imag(numpy.log(U_1)))

    return (U_1[inds],V_1[:,inds]) # = (U_d,V) s.t. U=V*U_d*V^\dagger




# ### $\phi_F$ from $U_d$
# Given $U_d$ by the function eigu, this function returns an array $\phi_F$.
# $$V e^{-i \phi_F} V^{\dagger} = V U_d V^{\dagger}$$ so, $$e^{-i \phi_F} = U_d$$ or $$ \phi_F = i \log(U_d) $$
# The following function returns and array vector $\phi_F$ for an operator as input.
# 

# In[50]:


def Phi_F(operator): 
    
    return (1j*np.log(eigu(operator)[0])).real  # eigu(Gr)[0] = exp(-i * phi_F).



## Ploting the graph

# In[68]:


f = open('plot_data.txt', 'w')

for i in range(1,N_eps):

    eps = i/(N_eps)
    
    f = open('plot_data'+Target_state+'.txt', 'a')
    
    X = str(eps)
    Y = Phi_F(Grover_reconstructed(eps))
            
    for j in range(2**N):
        f.write(X +','+ str(Y[j].real)+'\n')     

#print('Done')
#end = time.time()
#print('Time taken', (end-start)/60,'minutes')



#!/usr/bin/env python
# coding: utf-8

# In[131]:


import numpy as np


# In[132]:


# Number of qubits.
N = 4

L = N // 2 # Length of half cut number of qubits.

'''
Intiating a wave function with a lsit of size 2**N with all element as zeros.

''' 

Psi = [0]*(2**N)


# In[133]:


'''
Enter the non-zero coefficients of the wavefunction psi.

''' 

Psi[0] = 2/np.sqrt(13)
Psi[2] = np.sqrt(2/13)
Psi[4] = -1/np.sqrt(13)
Psi[14] = np.sqrt(6/13)


# Let $\Psi$ be the given normalized wave function. 
# Then $\rho_{AB} = <\Psi|\Psi>$ 
# is the density matrix of the subsystem $A$ and $B$, and 
# $\rho_{A} = Tr_{B}(\rho_{AB})$.
# The $(s,s^{'})$ element of the reduced density matrix $\rho_{A}$ is given by 
# $$<s|\rho_{A}|s^{'}> =\sum_{s^{''}} <ss^{''}|\rho_{AB}|s^{'}s^{''}>$$
# $$<s|\rho_{A}|s^{'}> =\sum_{s^{''}} <ss^{''}|\psi><\psi|s^{'}s^{''}>$$

# In[134]:


'''
It is not necessary to input a normalized wave function. 
It will be normalized later so that trace(rho) = 1.

'''

# Converting the list to a numpy matrix.
Psi = np.matrix(Psi).reshape(len(Psi),1) # Psi column matrix.

# Normalizing Psi.
Psi = Psi/np.linalg.norm(Psi)


# $$\psi_{s^{'}} = \Psi[ 2^{L} s^{'} : 2^{L} s^{'}+2^{L}-1] $$

# In[135]:


def psi(s):
    return Psi[(2**L)*s:(2**L)*s + 2**L]


# The elements of the matrix $\rho_{A}$ is given by
# $$<s|\rho_{A}|s^{'}> = \psi^{\dagger}_{s^{'}} \psi_{s}$$

# In[136]:


'''
    psi(s_p) is a row matrix/vector. psi(s) is a column matrix/vector.      
    Dimension of rhoA is N/2 x N/2. 
    The element <s|rhoA|sp> is given by psi_sp^\dagger * psi_s.
''' 

def rhoA(s,s_p): # <s|rho_A|s_p>
    


    # psi(s_p)^\dagger * psi(s) is the element of (s,s_p) of rho_AB.  
    return psi(s_p).getH() * psi(s)


# In[137]:


def rhoA_Matrix(N):
    
    M = np.zeros((N,N)) # 0 to N-1.
    
    '''
    rho is Hermitian, it is sufficient to calculate the elements above the
    diagonal. The the elements below the diagonal can be replace by the 
    complex cpnjugate of the elements above the diagonal.

    '''
    for i in range(N):
        for j in range(N):
            
            if i <= j : # Above the diagonal (i,j) i<j.
                
                M[i,j] = rhoA(i,j)[0,0]
                
            else: # Below the diagonal (i,j) i>j.
                
                M[i,j] = np.conjugate(M[j,i])
    return M


# In[138]:


'''
w is the diagonal of the diagonalized matrix rhoA.

'''
w, v = np.linalg.eig(rhoA_Matrix(N))


# In[140]:


DL = np.zeros(N) # Creating an array for log w with zeros.

'''
The following loop calculates S = - sum \lamba_i * log(\lambda_i).

'''

for i in range(len(w)):
    
    if abs(w[i]) < 1.e-8: # log of zero gives nan.
        
        pass # Leave the log(zero) element as zero.
    
    else:
        
        DL[i] = np.log(w[i])
        
# Entropy = -Tr(rho * log(rho)).        
print('Entropy S = ',-sum(w*DL))


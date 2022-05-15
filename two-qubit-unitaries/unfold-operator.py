import time
start = time.time()
import numpy as np
from two_level_unitary_unfold import Unfold
Target_state = '000000'
# First we note the length of N.
L = len(Target_state)

N = 2**L

# Then an identity matrix is created with the size 2**N with signs reversed.
U_w = - np.identity(2 ** L, dtype=complex)


Target_index = int(Target_state, 2)

## The position of the target state is set as 1.
U_w.itemset((Target_index, Target_index),1)

## We will first create a matrix with all elements 1. This is |psi><psi| = A(say).
A = np.ones((2**L, 2**L))

## U_s = 2\(2**N)2|psi><psi| - I
U_s = (2/(2**L))*A - np.identity(2**L)

## G is the Grover operator.
G = np.matmul(U_w, U_s)

print(len(Unfold(G)))

print("Time taken = ", time.time()-start, "seconds")

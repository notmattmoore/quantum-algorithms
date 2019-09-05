# Rough implementation of Simon's Algorithm using QuTiP
# Created by: Taylor Walenczyk
# Last Updated: 09.03.2019

from qutip import *
from itertools import *
import operator
import ualgebra as UA

# TODO Migrate utility functions to a seperate file/folder
# TODO Supress warnings about imaginary numbers

# An implementation of foldl (I think)
# In: func, the function to apply; xs, the list to fold
# Out: The result of folding the list with the specified function
def foldl(func, xs):
    ret = xs[0]
    for i in range(1,len(xs)):
        ret = func(ret, xs[i])
    return ret

# Converts integer to its k-bitstring representation
# In: x, an integer; k, the length of the bitstring
# Out: x as a k-bitstring
def int_to_bin(x, k):
    num = '{0:b}'.format(x)
    return '0'*(k-len(num))+num

def bin_to_int(x):
    ret = 0
    lenx = len(x)
    for i in range(lenx):
        ret += x[i]*(2**(lenx-1-i))
    return ret

# Generates list of binary strings, as specified
# In: k, desired length of representation; xs, list of integers to convert
# Out: A list of k-bitstrings in the order given
def gen_bin_list(k, xs):
    return [ int_to_bin(x,k) for x in xs ]

# Generates an oracle function for Simon's Algorithm
# In: k, order (?) of the group; D, the hidden subgroup D;
#     c, the common value for f(x) forall x in D
# Out: a function (array-form) that seperates cosets on the hidden subgroup
#     (f(a)=f(b) iff a-b in D)
def gen_oracle(k, D, c=0):
    bin_c = int_to_bin(c,k)
    return [ bin_c if int_to_bin(x, k) in D else int_to_bin(x, k) for x in range(2**k) ]
    # TODO improve ^ by minimizing "checks"

# Generates an oracle operator for Simon's Algorithm (Note: assumes structure
#     of oracle)
# In: k, the size of the group; f, the oracle function; mult, multiplier for the
#       if applicable
# Out: a (2**k)x(2**k) unitary operator embedding the oracle function
def gen_oracle_op(k, f, arity=0, mult=1):
    ret = [ [ 0 for _ in range(2**(k*arity)) ] for _ in range(2**(k*arity)) ]
    for x in range(len(f)):
        index = mult*(mult-1)
        for offset in range(len(f)):
            fx = bin_to_int(f[x])
            ket = tensor(   basis(2**k, mult-1),
                            basis(2**k, fx),
                            basis(2**k, (fx+offset) % 2**k) )
            # append just the array underlying the tensor
            ret[index+offset] = ket.full().astype(int).flatten().tolist()
    return ret

def gen_cong_op(k, A, cong):
    ret = [ basis(2**(2*k), i).full().astype(int).flatten().tolist() for i in range(2**(2*k)) ]
    for i,C in enumerate(UA.cong_classes(cong, A)):
        #out = int_to_bin(i,k)
        for x in C:
            xb = bin_to_int(x)
            for offset in range(2**k):
                ket = tensor( basis(2**k, xb), basis(2**k, (i+offset) % 2**k ))
                index = 2**k*xb+offset
                ret[index] = ket.full().astype(int).flatten().tolist()
    return ret

def gen_meet_op(n, A):
    meet_structure = [ [ meet(x,y) for y in A ] for x in A ]
    Meets = [
                Qobj(
                    inpt=gen_oracle_op(n, meet_structure[i], arity=3, mult=i+1),
                    dims=[ [ 2**n for _ in range(3) ] for _ in range(2) ])
                for i in range(len(meet_structure)) ]
    return foldl(operator.add, Meets)

# TODO Implement meet operator
# Simulates Simon's Algorithm
# In: k, order (?) of group k; P, an oracle operator encoding the homomorphism phi
# Out: TBD
def SemilatAlg(n,P,A):
    # Prepare the state psi
    G = [ basis(2**n, i) for i in range(2**n) ]
    psi = 2**(-n/2)*foldl(operator.add, G)
    print('~~~~~~~~ psi ~~~~~~~~')
    print(psi)

    # Prepare both registers
    regs = tensor(psi, basis(2**n,0))

    #print('~~~~~~~~ P ~~~~~~~~')
    #print(P)
    #print('~~~~~~~~ regs ~~~~~~~~')
    #print(regs)
    post_phi = P * regs
    #print('~~~~~~~~ post_phi ~~~~~~~~')
    #print(post_phi)

    # "discard" register 1
    #post_phi2 = post_phi.ptrace(1)
    #print('~~~~~~~~ post_phi2 ~~~~~~~~')
    #print(post_phi2)

    # augment state with |0>
    prep_meet = tensor( post_phi, basis(2**n, 0) )
    #prep_meet = tensor( post_phi2, psi, basis(2**n, 0) )

    # Apply meet operator
    M = gen_meet_op(n, A)
    # Meet acts on registers 1-3, so I prevents action on register 0
    #M_full = tensor(qeye(n),M)
    #print('~~~~~~~~ M ~~~~~~~~')
    #print(M)
    #print('~~~~~~~~ M_full ~~~~~~~~')
    #print(M_full)
    #print('~~~~~~~~ prep_meet ~~~~~~~~')
    #print(prep_meet)
    post_meet = M * prep_meet
    #print('~~~~~~~~ post_meet ~~~~~~~~')
    #print(post_meet)
    #print(post_meet.ptrace(0))

    #dm = post_meet * post_meet.dag()
    print('~~~~~~~~ dm final register partial ~~~~~~~~')
    print(dm.ptrace(2))



# ~~~ Testing ~~~

# Using the partial implementation of Simon's algorithm

def meet(x,y):
  return [x[i]*y[i] for i in range(len(x))]

M = UA.Operation(meet, 2, "meet")

n = 3
A = UA.FancySet( initial=[list(a) for a in list(product([0,1],repeat=n))] )
Theta = UA.rand_cong(A, [M], num_gen=1, Progress=False)

print('A:')
print(str(A))
print("Theta has classes:")
for C in UA.cong_classes(Theta, A):
  print(C)

Phi = Qobj( inpt=gen_cong_op(n, A, Theta), dims=[ [2**n,2**n] for _ in range(2) ] )

SemilatAlg(n,Phi,A)

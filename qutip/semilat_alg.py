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
    s = '0'*(k-len(num))+num
    return [ int(d) for d in s ]

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
            mult_ket = tensor([ basis(2, 1) if d == 1
                                else basis(2, 0) for d in int_to_bin(mult-1, k) ])
            fx_ket = tensor([ basis(2, 1) if d == 1 else basis(2, 0) for d in f[x] ])
            offset_ket = tensor([ basis(2, 1) if d == 1
                                else basis(2, 0) for d
                                in int_to_bin((fx+offset)%2**k, k) ])
            ket = tensor(   mult_ket, fx_ket, offset_ket )
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
                xb_ket = tensor([ basis(2, 1) if d == 1 else basis(2, 0) for d in x ])
                offset_ket = tensor([ basis(2,1)
                                        if d == 1
                                        else basis(2,0) for d in int_to_bin((i+offset)%2**k, k) ])
                ket = tensor( xb_ket, offset_ket)
                index = 2**k*xb+offset
                ret[index] = ket.full().astype(int).flatten().tolist()
    return ret

def gen_meet_op(n, A):
    meet_structure = [ [ meet(x,y) for y in A ] for x in A ]
    Meets = [
                Qobj(
                    inpt=gen_oracle_op(n, meet_structure[i], arity=3, mult=i+1),
                    dims=[[2]*3*n for _ in range(2) ])
                for i in range(len(meet_structure)) ]
    return foldl(operator.add, Meets)

# TODO Implement meet operator
# Simulates Simon's Algorithm
# In: k, order (?) of group k; P, an oracle operator encoding the homomorphism phi
# Out: TBD
def SemilatAlg(n,P,A):
    M = gen_meet_op(n, A)
    # Prepare the state psi
    # NOTE Tensoring a bunch of small basis vectors may seem silly
    # considering other options. However, this properly builds up
    # the tensor structure. This quality is extremely useful when working
    # with larger versions of built in operators (i.e. H^n)
    zn = tensor([ basis(2, 0) for _ in range(n) ])
    ht = hadamard_transform(n)
    psi = ht * zn

    # Prepare input registers for Phi
    regs = tensor(psi, zn)

    # Apply Phi
    post_phi = P * regs

    # "discard" register 1
    print('~'*10,'Experimenting with partial trace on output of Phi','~'*10)
    #phi_res = post_phi.ptrace([ 2*n - i for i in range(n,0,-1) ])

    #phi_res_transform = ht * phi_res
    #print('~'*10,'phi_res_transform','~'*10)
    #print(phi_res_transform)

    post_phi_transform = tensor(ht, identity([2 for _ in range(n)])) * post_phi
    post_transform_res = post_phi.ptrace([ 2*n - i for i in range(n,0,-1) ])
    print('~'*10,'post_transform_res','~'*10)
    print(post_transform_res)

    # Experimental stuff
    #outs = [ M * tensor( Qobj(inpt=row, dims=[[2]*n,[1]*n]), psi, zn ) for row in phi_res.full() ]
    #dm = foldl(operator.add, outs)
    #dm_last = dm.ptrace([ 3*n-i for i in range(n,0,-1) ])

    # Apply Hadamard gate
    #yld = dm_last * ht
    #print('~~~~~~~~ yld ~~~~~~~~')
    #print(yld)
    #print('~~~~~~~~ yld*zn ~~~~~~~~')
    #print(yld*zn)

    print('~'*10,'Experimenting with I_n tensor M','~'*10)

    # augment M
    M_aug = tensor( identity([2 for _ in range(n)]), M)

    # prep full register
    full_reg = tensor( post_phi, psi, zn ) # 2n x n x n

    # Apply M_aug
    full_res = M_aug * full_reg

    print('~'*10, 'full_res', '~'*10)
    print(full_res.data)

    # Restrict the space to area of interest
    interest = full_res.ptrace( [4*n - i for i in range(n,0,-1)] )
    print('~'*10, 'interest', '~'*10)
    print(interest)

    # Apply Hadamard gates
    interest_transform = ht * interest
    print('~'*10, 'interest_transform', '~'*10)
    print(interest_transform.data)

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

Phi = Qobj( inpt=gen_cong_op(n, A, Theta), dims=[[2]*2*n, [2]*2*n] )

SemilatAlg(n,Phi,A)

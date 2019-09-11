# Rough implementation of Simon's Algorithm using QuTiP
# Created by: Taylor Walenczyk
# Last Updated: 08/22/2019

from qutip import *
from itertools import *
from utils import *
import operator
import ualgebra as UA

# TODO Supress warnings about imaginary numbers

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
def gen_oracle_op(k, f, arity=2):
    ret = [ [ 0 for _ in range(2**(k*arity)) ] for _ in range(2**(k*arity)) ]
    for x in range(len(f)):
        for offset in range(len(f)):
            fx = bin_to_int(f[x])
            x_ket = tensor([ basis(2, d) for d in int_to_bin(x,k) ])
            offset_ket = tensor([ basis(2, d) for d in int_to_bin((fx+offset)%2**k, k) ])
            ket = tensor( x_ket, offset_ket )
            ret[x+offset] = ket.full().astype(int).flatten().tolist()
    return ret

# TODO Implement measurement phase
# Simulates Simon's Algorithm
# In: n, the digits for the binary representation of group elements;
#       U, an oracle operator
# Out: TBD
def SimonsAlg(n,U):
    # Prepare the state psi
    zn = tensor([ basis(2, 0) for _ in range(n) ])
    ht = hadamard_transform(n)
    psi = ht * zn

    # Prepare input registers for Phi
    regs = tensor(psi, zn)

    print(U)
    print(regs)

    # Apply the oracle
    full_state = U * regs

    # Reduce the space to the register of interest
    targ = full_state.ptrace([ i for i in range(n) ])


# ~~~ Testing ~~~

# Using the partial implementation of Simon's algorithm

n = 2
f = gen_oracle(n, [ list(product([[0],[1]], repeat=2)) ])
print(f)

# Programmatically generated operator
U = Qobj( inpt=gen_oracle_op(n,f), dims=[[2]*2*n, [2]*2*n])

#print("Expected operator")
#print(U1.data)
#print("Generated operator")
#print(Ua1.data)

SimonsAlg(n,U)

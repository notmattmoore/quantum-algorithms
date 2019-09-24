# Rough implementation of Simon's Algorithm using QuTiP

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
    #bin_c = int_to_bin(c,k)
    return [ c if x in D else x for x in range(2**k) ]
    # TODO improve ^ by minimizing "checks"

# Generates an oracle operator for Simon's Algorithm (Note: assumes structure
#     of oracle)
# In: k, the size of the group; f, the oracle function; mult, multiplier for the
#       if applicable
# Out: a (2**k)x(2**k) unitary operator embedding the oracle function
def gen_oracle_op(k, f, arity=2):
    ket_as_list = lambda ket: ket.full().astype(int).flatten().tolist()
    ret = []
    for x in range(len(f)):
        print('x',x)
        for offset in range(len(f)):
            print('offset',offset)
            #fx = bin_to_int(f[x])
            fx = f[x]
            print('fx', fx)
            x_ket = tensor([ basis(2, d) for d in int_to_bin(x,k) ])
            offset_ket = tensor([ basis(2, d) for d in int_to_bin((fx+offset)%2**k, k) ])
            ket = tensor( x_ket, offset_ket )
            ret.append(ket_as_list(ket))
    return ret

# TODO Implement measurement phase
# Simulates Simon's Algorithm
# In: n, the digits for the binary representation of group elements;
#       U, an oracle operator
# Out: TBD
def SimonsAlg(n,U):
    # Useful structures
    zn = tensor([ basis(2, 0) for _ in range(n) ])
    ht = hadamard_transform(n)

    # Prepare the state psi
    psi = ht * zn

    # Apply the oracle
    full_reg = tensor(psi, zn)
    post_oracle = U * full_reg

    # Reduce the space to the register of interest
    targ = post_oracle.ptrace([ i for i in range(n) ])

    # Return to the starting space
    targ = ht * targ

    return targ

# ~~~ Testing ~~~

# Using the partial implementation of Simon's algorithm

n = 2
f = gen_oracle(n, set([0,1]) )
#print(f)

# Programmatically generated operator
op = gen_oracle_op(n,f)
for row in op:
    print(row)
U = Qobj( inpt=op, dims=[[2]*2*n, [2]*2*n])

dm = SimonsAlg(n,U)
dm_to_hist(dm)

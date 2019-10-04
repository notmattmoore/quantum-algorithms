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
def gen_oracle(k, D):
    c = 0
    f = [ -1 for _ in range(2**k) ]
    for x in range(2**k):
        for y in range(2**k):
            if x^y in D:
                if x == y and f[x] == -1:
                    f[x] = c
                    c += 1
                else:
                    f[x] = f[y] = f[x] if x < y else f[y]
    return f


# TODO Implement measurement phase
# Simulates Simon's Algorithm
# In: n, the digits for the binary representation of group elements;
#       U, an oracle operator
# Out: TBD
def SimonsAlg(n,U):
    # registers and initial state
    reg1 = tensor([ basis(2, 0) for _ in range(n) ])
    reg2 = tensor([ basis(2, 0) for _ in range(n) ])
    psi0 = tensor(reg1, reg2)

    # the first set of gates
    In = tensor( [identity(2) for _ in range(n)] )
    HI = tensor( hadamard_transform(n), In )

    # do the Z2 inverse Fourier transform on the first register
    psi1 = HI * psi0

    # apply the oracle
    psi2 = U * psi1

    # do the Z2 Fourier transform (same as inverse for Z2) on the first register
    psi3 = HI * psi2

    # measure the first register and return the density matrix
    rho = psi3.ptrace(list(range(n)))
    return rho

# ~~~ Testing ~~~

def verify_oracle(f, U, n):
    cmp_kets = lambda x,y: ket_as_list(x) == ket_as_list(y)
    for x in range(len(f)):
        for y in range(len(f)):
            # We expect U|xy>=|x,y oplus f(x)>
            fx = f[x]
            yop = (y + fx) % 2**n
            expectation = tensor( int_to_ket(x,n), int_to_ket(yop,n) )

            # Find the reality
            inp = tensor( int_to_ket(x,n), int_to_ket(y,n) )
            reality = U * inp

            if not cmp_kets(expectation, reality):
                print('x y', x,y)
                print('expecation: val ket', yop, expectation)
                print('reality: val ket', '-1', reality)
                print('inp', inp)
                return
    print('Oracle verified')


# Using the partial implementation of Simon's algorithm

n = 4
f = gen_oracle(n, set([0,13]) )

# Programmatically generated operator
op = gen_oracle_op(n,f)
U = Qobj( inpt=op, dims=[[2]*2*n, [2]*2*n])
verify_oracle(f, U, n)

dm = SimonsAlg(n,U)
print(dm)
dm_to_hist(dm)

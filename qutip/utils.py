# Useful functions for implementing QuTiP algorithms

from qutip import *

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

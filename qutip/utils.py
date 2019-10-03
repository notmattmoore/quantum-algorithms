# Useful functions for implementing QuTiP algorithms

from qutip import *
import matplotlib.pyplot as plt

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

# Plots a density matrix as a histogram
# In:   dm, a density matrix
# Out:  nothing
def dm_to_hist(dm):
    names = [ str(x) for x in range(len(dm.diag())) ]
    values = dm.diag()
    plt.figure()
    plt.xlabel('Basis vector')
    plt.ylabel('Expectation value')
    plt.title('Probability Distribution')
    plt.bar(names, values)
    axes = plt.gca()
    axes.set_ylim([0.0,1.0])
    plt.show()

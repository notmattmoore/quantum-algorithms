# Useful summations
# Author: Taylor Walenczyk
# Last updated: 08.29.2019

import numpy as np
import operator
from .number_representation import arr_to_int

# In: a,b in ZZ_2^k as arrays
# Out: Dot product
def dot_prod(a, b): #{{{
    return np.dot(a,b)
# ----------------------------------------------------------------------------}}}
# solves the summation we get when evaluating
# Simon's algorithm on a semilattice
# In: group (list of k-bit strings), subgroup (list of k-bit strings in group)
def simonsSum(group, subgroup, k): #{{{
    mat = [[0 for _ in range(2**k)] for _ in range(2**k)] # Generates a 2^k dimensional matrix to match the group size (i.e. ZZ^2 has a 4-dimensional matrix)

    # a,b are members of the group
    for a in group:
        for b in group:

            # x,y are members of the hidden subgroup
            for x in subgroup:
                for y in subgroup:
                    xi = bin_to_int(x)
                    yi = bin_to_int(y)
                    if xi == yi:
                        continue
                    # print("({0},{1})".format(xi,yi))
                    mat[bin_to_int(a)][bin_to_int(b)] += (-1)**(abs(dot_prod(a,x)%2-dot_prod(b,y)%2))
        # end for
    #end for
    return mat
# ----------------------------------------------------------------------------}}}
# An implementation of foldl (I think)
# In: a function to apply func and a list to fold xs
# Out: The list folded according to the function
def foldl(func, xs): #{{{
    ret = xs[0]
    for i in range(1,len(xs)):
        ret = func(ret, xs[i])
    return ret
# ----------------------------------------------------------------------------}}}

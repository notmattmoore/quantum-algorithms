# Integer representation manipulations
# Author: Taylor Walenczyk
# Last updated: 08.29.2019

# splits a binary string into an array of integers
# In: x, a bitstring
# Out: x, in array form
def bin_to_arr(x): #{{{
    return [ int(d) for d in x ]
# ----------------------------------------------------------------------------}}}
# In: x, a bitstring stored in an array
# Out: the base 10 representation of x
def arr_to_int(x): #{{{
    ret = 0
    n = len(x)
    for i in range(n):
        ret += x[i]*2**(n-i-1)
    return ret
# ----------------------------------------------------------------------------}}}
# Ensures that the binary string in x has the proper number of digits in its representation
# In:   x, the binary number (in string form); n, the number of digits in the representation
# Out:  x represented by the proper number of binary digits in string form
def pad_bin(x,n): #{{{
    if len(x) > n:
        raise Exception('{0} is too long'.format(x))

    ret = x
    for _ in range(n-len(x)):
        ret = '0'+ret
    return ret
# ----------------------------------------------------------------------------}}}
# Higher order function making conversion to binary programmatically easier to read
# In:   x, the integer to convert to binary
def int_to_bin(x): #{{{
    return '{0:b}'.format(x)
# ----------------------------------------------------------------------------}}}
def int_to_arr(x,k):
    return bin_to_arr(pad_bin(int_to_bin(x),k))

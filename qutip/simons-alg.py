# Rough implementation of Simon's Algorithm using QuTiP
# Created by: Taylor Walenczyk
# Last Updated: 08/22/2019

from qutip import *
import operator

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
# In: k, the size of the group; f, the oracle function; c, the common value for
#     all elements of D, the hidden subgroup
# Out: a (2**k)x(2**k) unitary operator embedding the oracle function
def gen_oracle_op(k, f, c=0):
    ret = []
    for x in range(len(f)):
        for offset in range(len(f)):
            ket = tensor(   basis(2**k, x),
                            basis(2**k, (int(f[x], 2) + offset) % 2**k) )
            # append just the array underlying the tensor
            ret.append(ket.full().astype(int).flatten().tolist())
    return ret

# TODO Implement measurement phase
# Simulates Simon's Algorithm
# In: k, order (?) of group k; D, a hidden subgroup of G; U, an oracle operator
# Out: TBD
def SimonsAlg(k,U):
    # Prepare the state psi
    k = 2
    G = [ basis(2**k, i) for i in range(2**k) ]
    psi = 2**(-k/2)*foldl(operator.add, G)

    # Prepare both registers
    regs = tensor(psi, basis(2**k,0))

    # Apply U (after creating U)
    #print(U.data)
    #print(regs.data)

    #print(U.dims, U.shape, U.isherm)
    #print(regs.dims, regs.shape, regs.type)

    full_state = U * regs
    print(full_state.data)



# ~~~ Testing ~~~

# Using the partial implementation of Simon's algorithm

# Hard-coded operator
# U : G -> X s.t. U : |x,y> |-> |x,y+f(x)>
U1 = Qobj(
        inpt=[
            [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], # Identity on U|00,y>
            [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0], # ""
            [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0], # ""
            [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0], # ""
            [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0], # Identify on U|01,y>
            [0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0], # ""
            [0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0], # ""
            [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0], # ""
            [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0], # Permute bits on U|10,y> because f(x)=x
            [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0], # ""
            [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0], # ""
            [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0], # ""
            [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0], # Permute bits on U|11,y> because f(x)=x
            [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0], # ""
            [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1], # ""
            [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0]  # ""
        ],
        dims = [[4,4],[4,4]]
    )

D1 = gen_bin_list(2,[0,1])
f1 = gen_oracle(2,D1) # f(x) = x if x is not in D, else f(x) = 0 from 0-3

# Programmatically generated operator
Ua1 = Qobj(
            inpt=gen_oracle_op(2,f1)
        )

#print("Expected operator")
#print(U1.data)
#print("Generated operator")
#print(Ua1.data)

SimonsAlg(2,U1)

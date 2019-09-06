# Rough implementation of Semi-lattice Algorithm using QuTiP

from qutip import *
from itertools import *
from utils import *
import operator
import ualgebra as UA

# TODO Migrate utility functions to a seperate file/folder
# TODO Supress warnings about imaginary numbers


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

# Simulates the experimental Semilattice algorithm
# In: n, the number of digits needed to represent group elements in binary;
#       P, an oracle operator encoding the homomorphism phi
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

    # NOTE Overflow occurs for n > 3
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

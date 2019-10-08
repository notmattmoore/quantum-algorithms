# Rough implementation of Semi-lattice Algorithm using QuTiP

from qutip import *
from itertools import *
from utils import *
from random import *
import operator
import ualgebra as UA
import post_ops as PO

# TODO Supress warnings about imaginary numbers

def meet_set(S):
  # take the meet of a set of elements
  r = S[0]
  for s in S[1:]:
    r = meet(r, s)
  return r

def cong_to_oracle(n, Theta, A):
    f = [ x for x in range(2**n) ]
    for i,C in enumerate(UA.cong_classes(Theta,A)):
        for x in C:
            xb = bin_to_int(x)
            f[xb] = i
    return f

def gen_cong_op(n, A, cong):
    structure = [ [] for _ in range(2**(2*n)) ]

    # Encode the congruence classes
    for fx,C in enumerate(UA.cong_classes(cong, A)):
        for c in C:
            x = bin_to_int(c)
            for y in range(2**n):
                x_ket = int_to_ket(x, n)
                y_ket = int_to_ket((fx+y)%2**n, n)
                ket = tensor( x_ket, y_ket)
                for i,entry in enumerate(ket_as_list(ket)):
                    structure[i].append(entry)
    return structure

def phi_circuit(n,P):
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
    return post_phi

# Simulates the experimental Semilattice algorithm
# In: n, the number of digits needed to represent group elements in binary;
#       P, an oracle operator encoding the homomorphism phi
# Out: TBD
# TODO Determine where to go with this
def clone_alg(n,P):
    # Useful structures
    ht = hadamard_transform(n)

    # Apply the oracle
    post_oracle = phi_circuit(n, P)

    # Return to the starting space
    targ = ptrace_wrt_regs(post_oracle, [0], n)
    targ = ht*targ*ht

    return targ

# ~~~ Testing ~~~

# Using the partial implementation of Simon's algorithm

Ops = PO.named_clone("MP")
#Ops = PO.named_clone("DM")
#Ops = PO.named_clone("AP0")
#Ops = PO.named_clone("AP1")
#Ops = PO.named_clone("AP")
#Ops = PO.named_clone("MPT0inf")
#Ops = PO.named_clone("MPT1inf")
#Ops = PO.named_clone("MEETP")

print("Ops")
for op in Ops:
    print(op.name)
n = 3
Dn = [list(a) for a in list(product([0,1],repeat=n))]
#A, A_gens = UA.rand_subalg(Dn, Ops, Progress=False)
A = UA.FancySet( initial=Dn )
Theta, Theta_gens = UA.rand_cong(A, Ops, num_gen=randrange(n), Progress=False)

f = cong_to_oracle(n, Theta, A)
print(f)

print('A:')
print(str(A))
print("Theta has classes:")
for C in UA.cong_classes(Theta, A):
  print(C)

structure = gen_oracle_op(n, f)
Phi = Qobj( inpt=structure, dims=[[2]*2*n, [2]*2*n] )
#print(Phi)

dm = clone_alg(n,Phi)
print(dm)
print("DM trace", dm.tr())
dm_to_hist(dm)

# Rough implementation of Semi-lattice Algorithm using QuTiP

from qutip import *
from itertools import *
from utils import *
import operator
import ualgebra as UA

# TODO Supress warnings about imaginary numbers

def gen_cong_op(k, A, cong):
    # Build out full structure
    ret = [ basis(2**(2*k), i).full().astype(int).flatten().tolist() for i in range(2**(2*k)) ]

    # Encode the congruence classes
    for i,C in enumerate(UA.cong_classes(cong, A)):
        for x in C:
            xb = bin_to_int(x)
            for offset in range(2**k):
                xb_ket = tensor([ basis(2, d) for d in x ])
                offset_ket = tensor([ basis(2,d) for d in int_to_bin((i+offset)%2**k, k) ])
                ket = tensor( xb_ket, offset_ket)
                index = 2**k*xb+offset
                ret[index] = ket.full().astype(int).flatten().tolist()
    return ret

def gen_op(k, f, arity=0, mult=1):
    ret = [ [ 0 for _ in range(2**(k*arity)) ] for _ in range(2**(k*arity)) ]
    for x in range(len(f)):
        index = mult*(mult-1)
        for offset in range(len(f)):
            fx = bin_to_int(f[x])
            mult_ket = tensor([ basis(2, d) for d in int_to_bin(mult-1, k) ])
            fx_ket = tensor([ basis(2, d) for d in f[x] ])
            offset_ket = tensor([ basis(2, d) for d in int_to_bin((fx+offset)%2**k, k) ])
            ket = tensor( mult_ket, fx_ket, offset_ket )
            ret[index+offset] = ket.full().astype(int).flatten().tolist()
    return ret

def gen_meet_op(n, A):
    meet_structure = [ [ meet(x,y) for y in A ] for x in A ]
    Meets = [
                Qobj(
                    inpt=gen_op(n, meet_structure[i], arity=3, mult=i+1),
                    dims=[[2]*3*n for _ in range(2) ])
                for i in range(len(meet_structure)) ]
    return foldl(operator.add, Meets)

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
def SemilatAlg(n,P,A):
    # Set up
    M = gen_meet_op(n, A)
    zn = tensor([ basis(2, 0) for _ in range(n) ])
    ht = hadamard_transform(n)
    psi = ht * zn
    I_n = identity([2 for _ in range(n)])
    gen_last_reg = lambda a: [a*n - i for i in range(n,0,-1)] # a is arity of function

    post_phi = phi_circuit(n,P)
    ht_aug = tensor(ht, I_n)
    post_phi_transform = ht_aug * post_phi
    ppt_tr = post_phi_transform.ptrace(gen_last_reg(2))

    # augment M
    M_aug = tensor( identity([2 for _ in range(n)]), M)

    # prep full register
    full_reg = tensor( post_phi, psi, zn )

    # Apply M_aug
    full_res = M_aug * full_reg

    # Reduce to the register of interest
    targ = full_res.ptrace( [4*n - i for i in range(n,0,-1)] )

    # Apply Hadamard gates
    targ_transform = ht * targ
    print('~'*10, 'interest_transform', '~'*10)
    print(targ_transform.data)

# ~~~ Experiments ~~~
def exp1(n,P,A):
    print('~'*10,'Experiment 1','~'*10)
    print('Meet every element of A with the Hadamard transformed result of Phi')
    M = gen_meet_op(n, A)
    zn = tensor([ basis(2, 0) for _ in range(n) ])
    ht = hadamard_transform(n)
    psi = ht * zn
    I_n = identity([2 for _ in range(n)])
    gen_last_reg = lambda a: [a*n - i for i in range(n,0,-1)] # a is arity of function

    post_phi = phi_circuit(n,P)
    ht_aug = tensor(ht, I_n)
    post_phi_transform = ht_aug * post_phi
    ppt_tr = post_phi_transform.ptrace(gen_last_reg(2))

    print('~'*10,'post_phi','~'*10)
    print(post_phi)
    print('~'*10,'post_phi_transform','~'*10)
    print(post_phi_transform)
    print('~'*10,'post_phi_transform traced','~'*10)
    print(ppt_tr)


    # Meet the resultant with all element of t
    for t in A:
        print('~'*10,t,'~'*10)
        t_ket = tensor([ basis(2, i) for i in t ])

        # Control Calculation
        control_reg = tensor(t_ket,t_ket,zn)
        control_res = M * control_reg
        control_out = control_res.ptrace(gen_last_reg(3))
        print('t ^ t')
        print(control_out)

        # Experimental Calculation
        M_aug = tensor( identity([2 for _ in range(n)]), M )
        exp_reg = tensor(post_phi_transform, t_ket, zn) # 2n x n x n
        exp_res = M_aug * exp_reg
        exp_out = exp_res.ptrace(gen_last_reg(4))
        print('phi_out ^ t')
        print(exp_out)

def exp2(n,P,A):
    print('~'*10,'Experiment 2','~'*10)
    print('Meet every element of A with the untransformed result of Phi')
    M = gen_meet_op(n, A)
    zn = tensor([ basis(2, 0) for _ in range(n) ])
    ht = hadamard_transform(n)
    psi = ht * zn
    I_n = identity([2 for _ in range(n)])
    gen_last_reg = lambda a: [a*n - i for i in range(n,0,-1)] # a is arity of function

    post_phi = phi_circuit(n,P)

    # Meet the resultant with all element of t
    for t in A:
        print('~'*10,t,'~'*10)
        t_ket = tensor([ basis(2, i) for i in t ])

        # Control Calculation
        control_reg = tensor(t_ket,t_ket,zn)
        control_res = M * control_reg
        control_out = control_res.ptrace(gen_last_reg(3))
        print('t ^ t')
        print(control_out)

        # Experimental Calculation
        M_aug = tensor( identity([2 for _ in range(n)]), M )
        exp_reg = tensor(post_phi, t_ket, zn) # 2n x n x n
        exp_res = M_aug * exp_reg
        exp_out = exp_res.ptrace(gen_last_reg(4))
        print('phi_out ^ t')
        print(exp_out)

def exp3(n,P,A):
    print('~'*10,'Experiment 3','~'*10)
    print('Meet every element of A with the untransformed result of Phi')
    print('then Hadamard transform the result')
    M = gen_meet_op(n, A)
    zn = tensor([ basis(2, 0) for _ in range(n) ])
    ht = hadamard_transform(n)
    psi = ht * zn
    I_n = identity([2 for _ in range(n)])
    gen_last_reg = lambda a: [a*n - i for i in range(n,0,-1)] # a is arity of function

    post_phi = phi_circuit(n,P)

    # Meet the resultant with all element of t
    for t in A:
        print('~'*10,t,'~'*10)
        t_ket = tensor([ basis(2, i) for i in t ])

        # Control Calculation
        control_reg = tensor(t_ket,t_ket,zn)
        control_res = M * control_reg
        control_out = control_res.ptrace(gen_last_reg(3))
        print('t ^ t')
        print(control_out)

        # Experimental Calculation
        M_aug = tensor( identity([2 for _ in range(n)]), M )
        exp_reg = tensor(post_phi, t_ket, zn) # 2n x n x n
        exp_res = M_aug * exp_reg
        exp_out = exp_res.ptrace(gen_last_reg(4))
        exp_out_ht = ht * exp_out
        print('phi_out ^ t transformed')
        print(exp_out_ht)
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

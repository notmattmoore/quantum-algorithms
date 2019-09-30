# imports {{{1
from IPython import embed
from itertools import *
from random import *
from sys import *
import ualgebra as UA
import post_ops as PO
#----------------------------------------------------------------------------}}}1

def simon_sum(cong, a, b, verbose=False): # {{{
  dot = lambda x, y: sum([x[i]*y[i] for i in range(len(x))])
  exp = lambda x, y: dot(a,x) + dot(b,y)

  S = 0
  for [s, t] in cong:
    print("s t", s,t)
    S += (-1) ** (exp(s,t) % 2)

  if verbose:
    print("a b", a, b)
    print("-"*80)
    for [s, t] in cong:
      print("s, t:", s, t, "| a.s + b.t:", dot(a,s), "+", dot(b,t), "=", exp(s,t))
    print("+"*80)
    print("S", S)

  return S
#----------------------------------------------------------------------------}}}

# operations {{{
def meet(x,y):
  return [min(x[i], y[i]) for i in range(len(x))]
def join(x,y):
  return [max(x[i], y[i]) for i in range(len(x))]

def meet_set(S):
  # take the meet of a set of elements
  r = S[0]
  for s in S[1:]:
    r = meet(r, s)
  return r

def maj01(*args):
  args0_count = args.count(args[0])
  if args0_count >= len(args) - 1:
    return args[0]
  elif args0_count == 1 and args.count(args[1]) >= len(args) - 1:
    return args[1]
  return args[0]
def maj(*args):
  return list( map(maj01, *args) )

def extend_ops_cwise(Ops):
    new_Ops = []
    for op in Ops:
        temp_func = op.function
        new_op = UA.Operation(lambda *args : list( map(temp_func, *args) ),
                op.arity, "componenent-wise "+op.name)
        new_Ops.append(new_op)
    return new_Ops


Ops = []
#Ops.append(UA.Operation(meet, 2, "meet"))
#Ops.append(UA.Operation(join, 2, "join")) # comment out for semilattices
#Ops.append(UA.Operation(maj, 3, "maj"))
#Ops = PO.named_clone("MPT0inf")
Ops = extend_ops_cwise(PO.named_clone("DM"))

n = 4
Dn = [list(a) for a in product([0,1],repeat=n)]   # {0,1}^n

passes = True
while passes and count < 100:
  A, A_gens = UA.rand_subalg(Dn, Ops, Progress=False)
  #A = UA.FancySet( initial=Dn )
  #A_gens = []
  Theta, Theta_gens = UA.rand_cong(A, Ops, num_gen=randrange(n), Progress=False)
  min_preimages = [meet_set(C) for C in UA.cong_classes(Theta, A)]
  stdout.write(str(round(len(min_preimages) / len(A), 2)) + " ")
  stdout.flush()
  print("A", str(A))
  for a,b in product(A,repeat=2):
    S = simon_sum(Theta, a, b)
    if (S != 0 and a != b) or (S == 0 and a == b in min_preimages):
      passes = False
      break

print("A generators:")
for g in A_gens:
  print("  ", g)
print("Theta generators:")
for g in Theta_gens:
  print("  ", g)
print("Congruence classes:")
for C in UA.cong_classes(Theta, A):
  print("  ", C)
for a,b in product(A, repeat=2):
  S = simon_sum(Theta, a, b)
  if (S != 0 and a != b) or (S == 0 and a == b in min_preimages):
    S = simon_sum(Theta, a, b, verbose=True)
    print("a b S:", a, b, S)
    break

#print("Diagonal Entries")
#for a in A:
#    print('a', simon_sum(Theta, a, a))

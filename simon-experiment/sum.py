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
  return min(x,y)
def join(x,y):
  return max(x,y)

def meet_set(S):
  # take the meet of a set of elements
  r = [ min([ S[i][j] for i in range(len(S)) ]) for j in range(len(S[0])) ]
  return r

def NU(*args):
  args0_count = args.count(args[0])
  if args0_count >= len(args) - 1:
    return args[0]
  elif args0_count == 1 and args.count(args[1]) >= len(args) - 1:
    return args[1]
  return args[0]
def maj(*args):
  return list( map(NU, *args) )
#----------------------------------------------------------------------------}}}

Ops = []
#Ops.append(UA.Operation(meet, 2, "meet"))
#Ops.append(UA.Operation(join, 2, "join")) # comment out for semilattices
#Ops.append(UA.Operation(maj, 3, "maj"))

# PASSING
#Ops = extend_ops_cwise(PO.named_clone("DM"))
#Ops = extend_ops_cwise(PO.named_clone("MP"))
#Ops = extend_ops_cwise(PO.named_clone("AP0"))
#Ops = extend_ops_cwise(PO.named_clone("MPT0inf"))

# FAILING

# INDETERMINATE
#Ops = extend_ops_cwise(PO.named_clone("T"))
#Ops = extend_ops_cwise(PO.named_clone("P0"))
#Ops = extend_ops_cwise(PO.named_clone("P1"))
#Ops = extend_ops_cwise(PO.named_clone("P"))
#Ops = extend_ops_cwise(PO.named_clone("T0inf"))
#Ops = extend_ops_cwise(PO.named_clone("PT0inf"))
#Ops = extend_ops_cwise(PO.named_clone("T1inf"))
#Ops = extend_ops_cwise(PO.named_clone("PT1inf"))
#Ops = extend_ops_cwise(PO.named_clone("M"))
#Ops = extend_ops_cwise(PO.named_clone("MP0"))
#Ops = extend_ops_cwise(PO.named_clone("MP1"))
#Ops = extend_ops_cwise(PO.named_clone("MT0inf"))
#Ops = extend_ops_cwise(PO.named_clone("MPT02"))
#Ops = extend_ops_cwise(PO.named_clone("MT1inf"))
#Ops = extend_ops_cwise(PO.named_clone("MPT12"))
Ops = extend_ops_cwise(PO.named_clone("MPT1inf"))
#Ops = extend_ops_cwise(PO.named_clone("MEET"))
#Ops = extend_ops_cwise(PO.named_clone("MEETP0"))
#Ops = extend_ops_cwise(PO.named_clone("MEETP1"))
#Ops = extend_ops_cwise(PO.named_clone("MEETP"))
#Ops = extend_ops_cwise(PO.named_clone("JOIN"))
#Ops = extend_ops_cwise(PO.named_clone("JOINP0"))
#Ops = extend_ops_cwise(PO.named_clone("JOINP1"))
#Ops = extend_ops_cwise(PO.named_clone("JOINP"))
#Ops = extend_ops_cwise(PO.named_clone("D"))
#Ops = extend_ops_cwise(PO.named_clone("DP"))
#Ops = extend_ops_cwise(PO.named_clone("A"))
#Ops = extend_ops_cwise(PO.named_clone("AD"))
#Ops = extend_ops_cwise(PO.named_clone("AP1"))
#Ops = extend_ops_cwise(PO.named_clone("AP"))
#Ops = extend_ops_cwise(PO.named_clone("U"))
#Ops = extend_ops_cwise(PO.named_clone("UD"))
#Ops = extend_ops_cwise(PO.named_clone("UM"))
#Ops = extend_ops_cwise(PO.named_clone("UP0"))
#Ops = extend_ops_cwise(PO.named_clone("UP1"))
#Ops = extend_ops_cwise(PO.named_clone("F"))

for op in Ops:
    print("Name:", op.name)

n = 3
Dn = [list(a) for a in product([0,1],repeat=n)]   # {0,1}^n
min_preimages = []

passes = True
while passes:
  #A, A_gens = UA.rand_subalg(Dn, Ops, Progress=False)
  A, A_gens = UA.FancySet( initial=Dn ), []
  Theta, Theta_gens = UA.rand_cong(A, Ops, num_gen=randrange(n), Progress=False)
  min_preimages = UA.FancySet( initial=[meet_set(C) for C in UA.cong_classes(Theta, A)] )
  stdout.write(str(len(min_preimages)) + " ")
  stdout.flush()
  #for a,b in product(A,repeat=2):
  for a in A:
    b = a
    S = simon_sum(Theta, a, b)
    if (S != 0 and not (a == b in min_preimages)) or (S == 0 and a == b in min_preimages):
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
#for a,b in product(A, repeat=2):
for a in A:
  b = a
  S = simon_sum(Theta, a, b)
  if S != 0 and (a == b in min_preimages) or (a not in min_preimages):
    pass
  else:
    S = simon_sum(Theta, a, b, verbose=True)
    print("a b S:", a, b, S)
    break

#print("Diagonal Entries")
#for a in A:
#    print('a', simon_sum(Theta, a, a))

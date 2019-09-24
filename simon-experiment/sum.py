# imports {{{1
from IPython import embed
from itertools import *
from random import *
from sys import *
import ualgebra as UA
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

def meet(x,y):
  # meet as the bit-wise product
  return [x[i]*y[i] for i in range(len(x))]
def join(x,y):
  # join as the bit-wise truncated sum
  return [max(x[i],y[i]) for i in range(len(x))]

def meet_set(S):
  # take the meet of a set of elements
  r = S[0]
  for s in S[1:]:
    r = meet(r, s)
  return r

Ops = []
Ops.append(UA.Operation(meet, 2, "meet"))
Ops.append(UA.Operation(join, 2, "join")) # comment out for semilattices

n = 5
Dn = [list(a) for a in product([0,1],repeat=n)]   # {0,1}^n

passes = True
while passes:
  A, A_gens = UA.rand_subalg(Dn, Ops, Progress=False)
  Theta, Theta_gens = UA.rand_cong(A, Ops, num_gen=randrange(n), Progress=False)
  min_preimages = [meet_set(C) for C in UA.cong_classes(Theta, A)]
  stdout.write(str(round(len(min_preimages) / len(A), 2)) + " ")
  stdout.flush()
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

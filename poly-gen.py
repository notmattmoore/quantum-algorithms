# experimental framework to gather statistics for the probability of generating
# a lattice given a certain number of random generators

# imports {{{1
from itertools import *
from sys import *
import ualgebra as UA
#----------------------------------------------------------------------------}}}1

# experiment for distributive lattices {{{1
def meet(x,y):
  # meet as the bit-wise product
  return [x[i]*y[i] for i in range(len(x))]
def join(x,y):
  # join as the bit-wise truncated sum
  return [max(x[i],y[i]) for i in range(len(x))]

Ops = []
Ops.append(UA.Operation(meet, 2, "meet"))
Ops.append(UA.Operation(join, 2, "join"))

n = 10
A = UA.FancySet( initial=[list(a) for a in list(product([0,1],repeat=n))] ) # {0,1}^n

stats = dict()
for ng in range(n**2, 0, -1):
  stats[ng] = 0
  total=10**3
  for count in range(total):
    R, G = UA.rand_subalg(A, Ops, num_gen=ng, Progress=False)
    if len(R) == len(A):
      stats[ng] += 1
    stdout.write("\r" + str(ng) + ": " + str(round(count/total*100, 4)) + "%    ")
    if count % 10 == 0: stdout.flush()
  stdout.write("\n  " + "Generators = " + str(ng) + ", success = " + str(round(stats[ng]/total*100, 4)) + "%\n")
#----------------------------------------------------------------------------}}}1

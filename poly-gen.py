# experimental framework to gather statistics for the probability of generating
# a lattice given a certain number of random generators

# imports {{{1
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

n = 7
A = UA.FancySet( initial=[list(a) for a in list(product([0,1],repeat=n))] ) # {0,1}^n

stats = dict()
total=10**3
for count in range(total):
  R, G = UA.rand_subalg(A, Ops, num_gen=n**2, Progress=False)
  if len(R) not in stats:
    stats[len(R)] = 0
  stats[len(R)] += 1
  stdout.write(str(round(count/total*100, 4)) + "%    \r")
  if count % 10 == 0: stdout.flush()

print("----------")
for key in stats:
  print(key,stats[key]/(count+1)*100)
print("----------")
#----------------------------------------------------------------------------}}}1

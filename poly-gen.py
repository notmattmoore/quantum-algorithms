# experimental framework to gather statistics for the probability of generating
# a lattice given a certain number of random generators

# imports {{{1
from itertools import *
from sys import *
import ualgebra as UA
import post_ops as PO
#----------------------------------------------------------------------------}}}1

# experiment for distributive lattices {{{1
def meet(x,y):
  # meet as the bit-wise product
  return [x[i]*y[i] for i in range(len(x))]
def join(x,y):
  # join as the bit-wise truncated sum
  return [max(x[i],y[i]) for i in range(len(x))]

Names = [
"T"    , "P0"    , "P1"    , "P"    , "M"   , "MP0"   , "MP1"   , "MP"   ,
"MEET" , "MEETP0", "MEETP1", "MEETP", "JOIN", "JOINP0", "JOINP1", "JOINP",
"D"    , "DP"    , "DM"   , "A"   , "AD"    , "AP0"   , "AP1"   , "AP"   ,
"U"    , "UD"    , "UM"   , "UP0" , "UP1"   , "F"
] + [ "T0k for 2 <= k < inf" ] + [ "T0inf" ]   + \
    [ "PT0k for 2 <= k < inf" ] + [ "PT0inf" ]  + \
    [ "T1k for 2 <= k < inf" ] + [ "T1inf" ]   + \
    [ "PT1k for 2 <= k < inf" ] + [ "PT1inf" ]  + \
    [ "MT0k for 2 <= k < inf" ] + [ "MT0inf" ]  + \
    [ "MPT0k for 2 <= k < inf" ] + [ "MPT0inf" ] + \
    [ "MT1k for 2 <= k < inf" ] + [ "MT1inf" ]  + \
    [ "MPT1k for 2 <= k < inf" ] + [ "MPT1inf" ]


print("Possible clones:")
for name in Names:
    print("  ", name)
inp = input("Input a clone name: ")
print("Clone name", inp, "provided")
Ops = PO.named_clone(inp)
print("Corresponding operations:")
for op in Ops:
    print("  ", op.name) 

n = 5 
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

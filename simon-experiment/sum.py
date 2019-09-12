# imports {{{1
from itertools import *
from sys import *
from collections import Counter
from number_representation import *
import ualgebra as UA
import random
#----------------------------------------------------------------------------}}}1

def simon_sum(cong, a, b): # {{{
  exp = lambda x, y: sum([a[i]*x[i] + b[i]*y[i] for i in range(len(x))])
  S = 0
  for [s, t] in cong:
    S += (-1) ** (exp(s,t) % 2)
  return S
#----------------------------------------------------------------------------}}}

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
#----------------------------------------------------------------------------}}}1

n = 4
A = UA.FancySet( initial=[list(a) for a in list(product([0,1],repeat=n))] ) # {0,1}^n
Theta, G = UA.rand_cong(A, Ops, MaxNew=len(A), Progress=True)

print("Generators:")
for g in G:
  print("  ", g)
print("Congruence classes:")
for C in UA.cong_classes(Theta, A):
  print("  ", C)

print("\nSimons over first register")
print("a b S")
#for a,b in product(A,repeat=2):
#for a,b in combinations(A,2):
for a in A:
  S = simon_sum(Theta, a, a)
  print("  ",a, a, S)

print("\nSimons over second register")
print("t | size of congruence")
for i,C in enumerate(UA.cong_classes(Theta, A)):
    print("  ", i, len(C))

# experiment for semilattices {{{1
def meet(x,y):
  # meet as the bit-wise product
  return [x[i]*y[i] for i in range(len(x))]

Ops = []
Ops.append(UA.Operation(meet, 2, "meet"))
#----------------------------------------------------------------------------}}}1

# Sampling experiment
class Bag:
    '''Extends the Counter objec to better support Bag functions'''
    def __init__(self, elements):
        self.elements = Counter(elements)

    def stuff(self, elements):
        self.elements.update(elements)

    def grab(self):
        return random.choice(sorted(self.elements.elements()))

class ProbDist:
    '''Mimics a probability distribution'''
    def __init__(self, n, A, dm):
        self.n = n
        es = dict()
        for a in A:
            an = arr_to_int(a)
            es[an] = int(2**n * dm[an][an])
        self.dist = Bag(es)

    def sample(self):
        return int_to_arr(self.dist.grab(),self.n)

# create desnity matrix (sorta)
dm = [ [ 0 for _ in A ] for _ in A ]
for a in A:
    an = arr_to_int(a)
    dm[an][an] = 2**(-2*n) * simon_sum(Theta, a, a)
for row in dm:
    print(row)
pd = ProbDist(n,A,dm)
samples = 1000
print("\n{0} samples from pd".format(samples))
res = Counter()
for i in range(samples):
    res.update([ arr_to_int(pd.sample()) ])
print('\noutcome | ratio')
for outcome, count in res.items():
    oa = int_to_arr(outcome,n)
    print(oa, '|', count/samples)
#----------------------------------------------------------------------------}}}1

import ualgebra as UA
from itertools import *

def simon_sum(cong, a, b):
  exp = lambda x, y: sum([a[i]*x[i] - b[i]*y[i] for i in range(len(x))])
  S = 0
  for [s, t] in cong:
    S += (-1) ** (exp(s,t) % 2)
  return S

def meet(x,y):
  return [x[i]*y[i] for i in range(len(x))]

def join(x,y):
    return [x[i]|y[i] for i in range(len(x))]


M = UA.Operation(meet, 2, "meet")
J = UA.Operation(join, 2, "join")

n = 3
A = UA.FancySet( initial=[list(a) for a in list(product([0,1],repeat=n))] )
Theta = UA.rand_cong(A, [M,J], num_gen=1, Progress=False)
print("Theta")
print(Theta)

print("")

print("Theta has classes:")
for C in UA.cong_classes(Theta, A):
  print(C)

print("")

print("a b S | a Theta b?")
for a,b in combinations(A,2):
  S = simon_sum(Theta, a, b)
  print(a, b, S, "|", [a,b] in Theta)

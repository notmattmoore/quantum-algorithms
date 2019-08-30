# Group implementation using algebra library
# Author: Taylor Walenczyk
# Last updated: 08.30.2019

from utils.summations import dot_prod
from utils.number_representation import *
from utils.printing import *
comp_lib = __import__('computation-library')
import itertools
import math
import random
import sys

# TODO Refactor to use FancySet class
# TODO Refactor use of elements to accomodate set
class Group(): #{{{
    # Represents groups in ZZ_2**k where k is variable. The group operation at the
    # moment is component wise (i.e. bitwise) multiplication. Future iterations will
    # let this be definable.

    # Functions:
    # __init__(...): Provide the order of the group and the group operation
    #       along with your choice of elements
    # populate_M(self, elements): provide a list of elements with which to populate
    #       the group operation look up table (i.e. M)
    # m(self,x,y): returns the result of the group operation (i.e. M[x][y])

    # TODO Refactor to use FancySet
    def __init__(self, k, gop, e, elements=[]): #{{{
        self.elements = elements
        self.K = k
        self.E = e
        self.gop = gop
        self.Inv = dict()
        self.Op = dict()
        self.populate_Op(self.elements)
    # ========================================================================}}}
    def __len__(self): #{{{
        return len(self.elements)
    # ========================================================================}}}
    def __str__(self): #{{{
        present =   '|G| = {0}\n'.format(len(self))
        present +=  'G = {0}\n'.format(self.elements)
        present +=  'G.E = '+self.e()+'\n'
        present +=  'G.Op = '+pprint_dict(self.Op)
        return present
    # ========================================================================}}}
    # TODO Determine whether or not I need this function
    def populate_Op(self, es): #{{{
        for a in es:
            if a not in self.Op: # Initialize if necessary
                self.Op[a] = dict()
            for b in self.elements:
                if b not in self.Op:
                    self.Op[b] = dict()
                if self.gop == None: # Who knows what to do then?
                    continue
                res = self.gop(a,b)
                add = pad_bin(int_to_bin(res), self.K)
                self.Op[a][b] = add
                self.Op[b][a] = add
                # looking forward to prevent issues
                if add not in self.Op:
                    self.Op[add] = dict()
    # ========================================================================}}}
    def op(self,x,y): #{{{
        return self.Op[x][y]
    # ========================================================================}}}
    def inv(self,x): #{{{
        if x in self.Inv:
            return self.Inv[x]
        for y in self.elements:
            if self.op(x,y) == self.e():
                self.Inv[x] = y
                return y
        return None # if x does not have an inverse
    # ========================================================================}}}
    def e(self): #{{{
        return self.E
    # ========================================================================}}}
    def is_group(self): #{{{
        # idempotence
        for x in self.elements:
            if self.op(x,self.e()) != x:
                return False

        # associativity
        for tup in itertools.product(self.elements,repeat=3):
            x = tup[0]
            y = tup[1]
            z = tup[2]

            try:
                if self.op(self.op(x,y),z) != self.op(x,self.op(y,z)):
                    print('NON-ASSOCIATIVE: x={0}, y={1}, z={2}'.format(x,y,z))
                    return False
            except:
                print('Invalid indexing: x={0}, y={1}, z={2}'.format(x,y,z))
                sys.exit()
        return True
    # ========================================================================}}}
    def is_abelian(self): #{{{
        # commutativity
        for tup in itertools.product(self.elements, repeat=2):
            x = tup[0]
            y = tup[1]

            try:
                if self.Op(x,y) == self.Op(y,x):
                    return False
            except:
                print('Invalid indexing: x={0}, y={1}'.format(x,y))
                sys.exit()
        return False
    # ========================================================================}}}
# ----------------------------------------------------------------------------}}}

# Creates a passable group function
# In:   G, a group
# Out:  a callable function that performs an embedded group's operation
def grp_op(G):
    def op(x,y):
        return G.op(x,y)
    return op

# Generates a proper group from generating elements
# In:   G, a group; gens, the generating elements
# Out:  a FancySet of elements in the subgroup
def gen_subgroup(G, gens):
    cur_set = comp_lib.FancySet(initial=gens, addl=[ 'generator' for _ in gens ])
    new_set = comp_lib.FancySet(initial=[G.e()], addl='identity element') # ensures e exists
    # Add inverse elements
    for x in cur_set:
        xi = G.inv(x)
        if xi and xi not in cur_set: # prevents overwriting addl
            new_set.add(xi, addl='Inverse of {0}'.format(x))

    while len(new_set) != 0: # while the closure grows
        cur_set.update(new_set)
        new_set = comp_lib.FancySet()
        for ps in itertools.product(cur_set, repeat=2):
            res = G.op(ps[0],ps[1])
            if res not in cur_set:
                note = 'Generated by G.op({0},{1})'.format(ps[0],ps[1])
                print(res+' '+note)
                new_set.add(res, addl='Generated by G.op({0},{1})'.format(ps[0],ps[1]))
        # MAY NEED TO WRAP THIS IN ITS OWN CLOSURE
        # Add inverse elements
        for x in new_set: # Presumes inverse are present already for old set
            xi = G.inv(x)
            if xi and xi not in new_set: # prevents overwriting addl
                new_set.add(xi, addl='Inverse of {0}'.format(x))
        # Add elements to ensure assoc
        new_union = cur_set.union(new_set)
        for tup in itertools.product(new_union, repeat=3):
            x = tup[0]
            y = tup[1]
            z = tup[2]

            a = G.op(x,y)
            b = G.op(y,z)
            c = G.op(a,z)

            if a not in cur_set:
                new_set.add(a)
            if b not in cur_set:
                new_set.add(b)
            if c not in cur_set:
                new_set.add(c)
        # Add inverse elements
        for x in new_set:
            xi = G.inv(x)
            if xi and xi not in new_set: # prevents overwriting addl
                new_set.add(xi, addl='Inverse of {0}'.format(x))
    return cur_set

# Creates a random subgroup
# In:   num_gens, the number of generators for the subgroup (-1 if a random amount
#       is desired); G, a group object; Ops, list of relevant operators
# Out:  whoknows
def rand_subgroup(G, num_gens=-1):
    if num_gens == -1:
        num_gens = random.randrange(len(G))
    gens = comp_lib.FancySet()
    for _ in range(num_gens):
        gens.add(random.choice(G.elements), addl='generator')
    S = gen_subgroup(G, gens)
    return S

g = Group(
        k=3,
        gop=lambda x,y: int(x,2) & int(y,2),
        e='1'*3,
        elements=[ ''.join(tup) for tup in itertools.product( ['0','1'], repeat=3 ) ]
        )

print(str(g))

def tight_grp_op(x,y):
    return g.op(x,y)

gop = comp_lib.Operation(tight_grp_op, 2, 'componentwise mult')
sg = rand_subgroup(G=g)
print(str(sg))

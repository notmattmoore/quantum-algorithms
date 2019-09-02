# Version: 2019-09-02
# A library for doing computations on universal algebras

# imports {{{1
from copy import copy, deepcopy
from IPython import embed
from multiprocessing import Pool
from itertools import *
from random import choice, randrange
from sys import *
import cProfile, pickle, random, shelve, string
#---------------------------------------------------------------------------}}}1

# data structures, etc
class FancySet: # {{{1
  # A class that extends the python set datatype to arbitrary (maybe unhashable)
  # elements. We use a dictionary with keys being the string representation. The
  # dictionary values will be pairs with the first coordinate being the stored
  # item and the second coordinate being an "additional information" field.
  # Aside from consolidating things, the main benefit is that the statement
  # 'x in FS' has average O(1) time complexity.

  # FS.elements: the dictionary

  # FS.__init__(...) allows you to specify the initial set and initial addl info.

  # FS.addl(e) returns the additional info for element e.
  # FS.add(e,addl=None) adds element e to FS. You can provide your addl
  #   information.
  # FS.remove(e) removes element e from FS, raising KeyError if e isn't found.
  # FS.union(T) returns a new FS that is the union of the two.
  # FS.update(T) adds elements of T to FS.

  # FS.__contains__(...) is called for 'x in FS'. Works based on the keys, so
  #   O(1) average time complexity.
  # FS.__leq__(...) is called for S <= T. Looks only at the keys.
  # FS.__eq__(...) is called for S == T. Looks only at the keys.

  # FS.__iter__(...) is called for 'iter(FS)' or implict things involving loops,
  #   etc. This gives the actual element values, not the strings which index
  #   them.
  # FS.__getitem__(...) is called for FS[...]. It returns the element.
  # FS.__len__(...) is called for 'len(FS)'.
  # FS.__str__(...) is called for 'str(FS)'.

  def __init__(self, initial=[], addl=[]): # {{{
    self.elements = dict()
    if addl == []:
      addl = [ None for _ in initial ]
    for i,e in enumerate(initial):
      self.add(e, addl[i])
  #--------------------------------------------------------------------------}}}

  def addl(self, e): # {{{
    return self.elements[str(e)][1]
  #--------------------------------------------------------------------------}}}

  def add(self, e, addl=None):   # {{{
    self.elements[str(e)] = (e, addl)
  #--------------------------------------------------------------------------}}}
  def remove(self, e):   # {{{
    del self.elements[str(e)]
  #--------------------------------------------------------------------------}}}
  def union(self, T):   # {{{
    S = deepcopy(self)
    S.update(T)
    return S
  #--------------------------------------------------------------------------}}}
  def update(self, T):   # {{{
    self.elements.update(T.elements)
  #--------------------------------------------------------------------------}}}

  def __contains__(self, e):  # {{{
    return str(e) in self.elements
  #--------------------------------------------------------------------------}}}
  def __leq__(self, T):  # {{{
    for k in self.elements:
      if k not in T:
        return False
    return True
  #--------------------------------------------------------------------------}}}
  def __eq__(self, T):  # {{{
    return len(self) == len(T) and self <= T
  #--------------------------------------------------------------------------}}}

  def __iter__(self):  # {{{
    for k in self.elements:
      yield self.elements[k][0]
  #--------------------------------------------------------------------------}}}
  def __len__(self):  # {{{
    return len(self.elements)
  #--------------------------------------------------------------------------}}}
  def __str__(self):  # {{{
    ret = ""
    for e in self.elements:
      ret += e + ",\n"
    return ret[:-2]
  #--------------------------------------------------------------------------}}}
#----------------------------------------------------------------------------}}}1
class ishelf(object):  # {{{1
# Class to store an array backed by a file using shelve. There is a key
# self.store['length'] which stores the length of an "array" (dictionary,
# actually) self.store['i']. Other keys can be stored to as well.
#
# IS.__init__(...) open the given filename with optional mode: r (read),
#   w (read+write existing), n (new/replace),
#   c (read+write, create if needed, default).
# IS.__del__(key) is called for 'del IS[index]'. We just close the shelve.
#
# IS.len(full=False) returns the length of the array unless full=True.
# IS.get(key,default=None) returns IS[key] if it exists, and default otherwise.
# IS.append(value, key=None) append with optional key, update length if no key
#   is provided
# IS.contains(value, full=False) T/F if value is in the *array*. Optionally
# IS.has_key(key) returns whether the given key is in the dictionary
# IS.keys() returns the keys of the dictionary.
# IS.iter(full=False) iterate over the values indexed values of A, or optionally
#   over over all values.
#
# IS.__len__() is called for 'len(IS)'. Return the length of the *array*.
# IS.__getitem__(key) is called to return 'IS[index]'
# IS.__setitem__(key, value) is called to assign 'IS[index]=value'. Does not
#   update length.
# IS.__delitem__(key) is called for 'del IS[index]'. We only support deletion of
#   string keys.
# IS.__contains__(value) is called for 'x in IS'. Only works for indexed values.
# IS.__iter__() is called for 'iter(IS)'. Iterates only over indexed values.
# IS.__str__() is called for 'str(IS)' and for 'print(IS)'. Only shows indexed
#   values.
#
# IS.close() close the file
# IS.sync() sync the shelf to the file

  store = None

  def __init__(self, filename, mode="c"): # {{{
    self.store = shelve.open(filename,flag=mode)
  #==========================================================================}}}
  def __del__(self):  # {{{
    self.close()
  #==========================================================================}}}

  def len(self, full=False): # {{{
    if not full:
      length = self.store.get('__length', 0)
    else:
      length = len(self.store)
    return length
  #==========================================================================}}}
  def get(self, key, default=None): # {{{
    if type(key) == int and key < self.len():
      return self[key]
    elif type(key) == str and key in self.keys():
      return self[key]
    return default
  #==========================================================================}}}
  def append(self, value, key=None): # {{{
    if key == None:
      key = self.len()
      self.store['__length'] = key + 1
    self[key] = value
  #==========================================================================}}}
  def extend(self, values): # {{{
    for v in values: self.append(v)
  #==========================================================================}}}
  def contains(self, value, full=False): # {{{
    if not full:
      keys = range(len(self))
    else:
      keys = self.store.keys()

    for i in keys:
      if self[i] == value: return True
    return False
  #==========================================================================}}}
  def has_key(self, key): # {{{
    return key in self.keys()
  #==========================================================================}}}
  def keys(self): # {{{
    return sorted( self.store.keys() )
  #==========================================================================}}}
  def iter(self, full=False): # {{{
    if not full:
      for i in range(self.len()): yield self[i]
    else:
      for key in self.keys(): yield (key, self.store[key])
  #==========================================================================}}}

  def __len__(self): # {{{
    return self.len()
  #==========================================================================}}}
  def __getitem__(self, key): # {{{
    if type(key) == int and not ( -len(self) <= key < len(self) ):
      raise IndexError("ishelf index out of range")

    if key < 0: key += len(self)

    return self.store.get(str(key))
  #==========================================================================}}}
  def __setitem__(self, key, value): # {{{
    self.store[str(key)] = value
  #==========================================================================}}}
  def __delitem__(self, key): # {{{
      if type(key) == str:
        del self.store[key]
      else:
        raise TypeError("ishelf does not support deletion of non-string keys.")
  #==========================================================================}}}
  def __contains__(self, value): # {{{
    return self.contains(value)
  #==========================================================================}}}
  def __iter__(self): # {{{
    return self.iter()
  #==========================================================================}}}
  def __str__(self): # {{{
    s = "["
    for i in range(len(self)-1): s += str(self[i]) + ", "
    return s + str(self[-1]) + "]"
  #==========================================================================}}}

  def close(self): # {{{
    self.store.close()
  #==========================================================================}}}
  def sync(self): # {{{
    self.store.sync()
  #==========================================================================}}}
#---------------------------------------------------------------------------}}}1

def save(S,filename,mode="w"):  # {{{
  f = open(filename+".pickle", mode)
  pickle.dump(S, f)
  f.close()
#----------------------------------------------------------------------------}}}
def load_one(filename):  # {{{
  f = open(str(filename)+".pickle","r")
  S = pickle.load(f)
  f.close()
  return S
#----------------------------------------------------------------------------}}}
def load_all(filename):  # {{{
  f = open(filename+".pickle", "r")
  while True:
    try:
      P = pickle.load(f)
      yield P
    except EOFError: break
  f.close()
#----------------------------------------------------------------------------}}}

# algebraic functions
class Operation: # {{{1
  # an operation of an algebra

  # O.function -- function itself
  # O.arity    -- arity of the function
  # O.name     -- the function symbol

  # We can call O(input), where input is a list of vectors in the domain of the
  # function (to some power). For example, O([[1,2],[3,4],[5,6]]) represents input
  # to a 3-ary function of 2-ary vectors: [O(1,3,5), O(2,4,6)]. This is hackish,
  # but is a workaround to pool.map not supporting multiple iterables.

  def __init__(self, function, arity, name):  # {{{
    self.function = function
    self.arity = arity
    self.name = name
  #--------------------------------------------------------------------------}}}
  def __call__(self, args):  # {{{
    return list( map(self.function, *args) )
  #--------------------------------------------------------------------------}}}
  def pprint(self, *args):  # {{{
    result = self(args)
    for row in range(len(args[0])):
      if row == 0:
        R = self.name
      else:
        R += "\n" + " "*len(self.name)
      R += "( "
      R += ", ".join( [ str(args[col][row]) for col in range(len(args)) ] )
      R += " ) "
      if row == 0:
        R += "="
      else:
        R+= " "
      R += " " + str(result[row])
    return R
#----------------------------------------------------------------------------}}}
#----------------------------------------------------------------------------}}}1

def is_reln(R, Ops, Skip = [], Progress=True, return_on_fail=True): # {{{
  cores = 4
  pool = Pool(cores)
  op_total = len(Ops)
  for op_count, op in enumerate(Ops):
    if op.name in Skip:
      continue
    args_total = len(R)**op.arity
    args_all = product(R, repeat=op.arity)
    # We have to seek the arguments that gave a certain result. In case we don't
    # return_on_fail, we should keep track of the previous index. Not a good
    # soln...
    args_all_seek = product(R, repeat=op.arity)
    prev_args_count = -1
    for args_count, result in enumerate(pool.imap(op, args_all, chunksize=1000)):
      if result not in R:
        stdout.write("\n\nFound " + str(result) + ":\n")
        args = next(islice(args_all_seek, args_count-prev_args_count-1, None))
        prev_args_count = args_count
        stdout.write(op.pprint(args) + "\nwhich is not in the relation.\n")
        if return_on_fail:
          return False
      if Progress and args_count % 100000 == 0:
        stdout.write( "\roperation " + op.name \
            + " " + str(op_count+1) + " / " + str(op_total) \
            + ", argument " + str(args_count) + " / " + str(args_total) \
            + " ~ " + str( round( args_count/args_total*100, 4 ) ) + "%  " )
        stdout.flush()
  if Progress: stdout.write("  done.\n")
  return True
#----------------------------------------------------------------------------}}}

def powerset_as_indicators(size): # {{{
  I = [0]*size
  T = [1]*size
  while I != T:
    yield copy(I)
    index = 0
    while I[index] == 1:
      I[index] = 0
      index += 1
    I[index] = 1
  yield copy(I)
#----------------------------------------------------------------------------}}}
def single_closure(G_old, G_new, Ops, Progress=True, Search=None):  # {{{
  # G_old is a set of elements. G_new has been computed by taking G_old and
  # applying single functions from Ops to it. It should be disjoint from G_old.
  # We return G_newer, which is the result of applying functions from Ops to
  # inputs from G_old and G_new, with at least 1 element of G_new as input. It
  # will be disjoint from G_old+G_new.

  cores = 4
  pool = Pool(cores)
  G_newer = FancySet()
  op_total = len(Ops)
  for op_count, op in enumerate(Ops):
    args_total = (len(G_old)+len(G_new))**op.arity - len(G_old)**op.arity
    args_count = 0
    # We don't want to evaluate op(all G_old), so we do all combinations of
    # elements from G_old and G_new, with at least 1 G_new always appearing. We
    # do this by quantifying over all variable positions where a new element
    # will appear.
    all_old = [0]*op.arity
    for vars_old_new in powerset_as_indicators(op.arity):
      if vars_old_new == all_old:
        continue
      args_all = product( *[ [G_old, G_new][var] for var in vars_old_new ] )
      # we have to seek the arguments that gave us result... not a good soln...
      args_all_seek = product( *[ [G_old, G_new][var] for var in vars_old_new ] )
      prev_args_index = -1
      for args_index, result in enumerate(pool.imap(op, args_all, chunksize=1000)):
        # if result is something new
        if not ( result in G_old or result in G_new or result in G_newer ):
          args = next(islice(args_all_seek, args_index-prev_args_index-1, None))
          prev_args_index = args_index
          G_newer.add(result, op.pprint(*args))
          if Search != None and Search(result):
            stdout.write( "\nFound:\n" + op.pprint(*args) + "\n" )
        args_count += 1
        if Progress and args_count % 10000 == 0:
          stdout.write( "\roperation " + op.name \
              + " " + str(op_count+1) + " / " + str(op_total) \
              + ", argument " + str(args_count) + " / " + str(args_total) \
              + " ~ " + str( round( args_count/args_total*100, 4 ) ) + "%" \
              + ", new elements: " + str(len(G_newer)) + " "*2 )
          stdout.flush()
  if Progress:
    stdout.write("  done.\n")
  return G_newer
#----------------------------------------------------------------------------}}}
def subalg_gen(Generators, Ops, Progress=True, Search=None, ExtraClosure=None, SavePartial=None, MaxLevels=-1):  # {{{
  G_old = FancySet()
  G_new = FancySet(initial=Generators)
  closure_level = 0
  if Progress:
    stdout.write("Generators:\n")
    for g in Generators:
      stdout.write("  " + str(g) + "\n")
  while len(G_new) != 0:
    closure_level += 1
    if Progress:
      stdout.write( "Closure level " + str(closure_level) + ", " )
      stdout.write( "at " + str(len(G_old) + len(G_new)) + " elements:\n" )
      stdout.flush()
    G_newer = single_closure(G_old, G_new, Ops, Progress=Progress, Search=Search)
    if ExtraClosure != None:
      G_extra = ExtraClosure(G_old, G_new, Ops, Search=Search)
      if Progress:
        stdout.write( "Extra closure found " + str(len(G_extra)) )
        stdout.write( " new elements.\n" )
      G_newer.update(G_extra)
    G_old.update(G_new)
    G_new = G_newer
    if SavePartial != None:
      save(G_old.union(G_newer), SavePartial + "_partial" + str(closure_level))
    if closure_level == MaxLevels:
      G_old.update(G_new)
      return G_old

  return G_old
#----------------------------------------------------------------------------}}}
def subalg_gen_layers(Generators, Ops, Progress=False):  # {{{
  G_old = FancySet()
  G_new = FancySet(initial=Generators)
  closure_level = 0
  if Progress:
    stdout.write("Generators:\n")
    for g in Generators:
      stdout.write("  " + str(g) + "\n")
  while len(G_new) != 0:
    yield G_new
    closure_level += 1
    if Progress:
      stdout.write( "Closure level " + str(closure_level) + ", " )
      stdout.write( "at " + str(len(G_old) + len(G_new)) + " elements:\n" )
      stdout.flush()
    G_newer = single_closure(G_old, G_new, Ops, Progress=Progress)
    G_old.update(G_new)
    G_new = G_newer
#----------------------------------------------------------------------------}}}

def transitive_closure_layer(C, C_new, A, Search=False):  # {{{
  # C and C_new are sets of elements, with C_new disjoint from C. A is the
  # underlying algebra. C is assumed to be transitively closed. Returns a set
  # C_newer of elements that are not in C u C_newer but are in the transitive
  # closure. 

  C_newer = FancySet()

  # C is transitive, so we don't need to check C x C. We do need to check
  # combinations of C x C_new and C_new x C_new. Note that the transitive
  # via C x C_new produces the same as via C_new x C. The outer loop controls
  # which of these combinations we're looking at.
  for [comb1, comb2] in [ (0,1), (1,1) ]:
    C_comb1 = [C, C_new][comb1]
    C_comb2 = [C, C_new][comb2]
    for [a,b], c in product(C_comb1, A):
      if [b,c] in C_comb2 and not ([a,c] in C or [a,c] in C_new or [a,c] in C_newer):
        why = "TransitiveClosure( " + str([a,b]) + ", " + str([b,c]) + " ) = " + str([a,c])
        C_newer.add([a,c], why)
        if Search != None and Search([a,c]):
          stdout.write( "\nFound:\n" + why + "\n" )
  return C_newer
#----------------------------------------------------------------------------}}}
def cong_gen(Generators, Ops, Progress=True, Search=None, SavePartial=None, MaxLevels=-1):  # {{{

  A = FancySet()
  GeneratorsCong = FancySet( initial=Generators )

  # Make sure that Generators contains the diagonal and is symmetric
  for [a,b] in Generators:
    GeneratorsCong.add([b,a])
    GeneratorsCong.add([a,a])
    GeneratorsCong.add([b,b])
    A.add(a)
    A.add(b)

  # subalg_gen expects the ExtraClosure function to take certain arguments, so
  # we make the function we want to send it look like that
  def transitive_closure_wrapper(C, C_new, Ops, A=A, Search=Search):
    return transitive_closure_layer(C, C_new, A, Search=Search)

  return subalg_gen(GeneratorsCong, Ops, Progress=Progress, Search=Search, \
      ExtraClosure=transitive_closure_wrapper, SavePartial=SavePartial, \
      MaxLevels=MaxLevels)
#----------------------------------------------------------------------------}}}
def rand_cong(A, Ops, num_gen=-1, Progress=True): # {{{
  if num_gen == -1:
    num_gen = randrange(1,len(A)+1)

  G = [ [a,a] for a in A ]
  A_list = list(A)
  for _ in range(num_gen):
    a = choice(A_list)
    b = choice(A_list)
    G.append([a,b])
  return cong_gen(G, Ops, Progress=Progress)
#----------------------------------------------------------------------------}}}
def cong_classes(C, A): # {{{
  # output the congruence classes of C
  found = FancySet()
  classes = []
  for a in A:
    if a not in found:
      classes.append([])
      for b in A:
        if [a,b] in C:
          found.add(b)
          classes[-1].append(b)
  return classes
#----------------------------------------------------------------------------}}}


#Example:
#def m(x,y):
#  z = [-1]*len(x)
#  for i in range(len(x)):
#    z[i] = x[i]*y[i]
#  return z
#
#M = Operation(m, 2, "m")
#
#A = FancySet( initial=[list(a) for a in list(product([0,1],repeat=3))] )
#G = [ [[1,1,1], [1,0,1]] ] + [ [a]*2 for a in A ]
#C = cong_gen(G, [M])
#D = rand_cong(A, [M], 1)
#
#found = []
#for a in A:
#  if a in found:
#    continue
#  for b in A:
#    if b in found:
#      continue
#    if [a,b] in D:
#      found.append(b)
#      stdout.write(str(b) + " ")
#  stdout.write("\n")

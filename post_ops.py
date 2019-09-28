"""
Library providing access to operational bases for each of the clones in Post's
lattice. See http://en.wikipedia.org/wiki/Post's_lattice.
"""

# Name            Basis                                                      {{{
# ------------------------------------------------------------------------------
# T               join, neg
# P0              join, +
# P1              meet, implies
# P               x ? y : z (ternary conditional, q)
# T0k, k >= 2     th_k^k+1, not_implies
# T0inf           not_implies
# PT0k, k >= 2    th_k^k+1, x meet (y implies z)
# PT0inf          x meet (y implies z)
# T1k, k >= 2     th_2^k+1, implies
# T1inf           implies
# PT1k, k >= 2    th_2^k+1, x join (y + z)
# PT1inf          x join (y + z)
# M               meet, join, 0, 1
# MP0             meet, join, 0
# MP1             meet, join, 1
# MP              meet, join
# MT0k, k >= 2    th_k^k+1, 0
# MT0inf          x meet (y join z), 0
# MPT0k, k >= 3   th_k^k+1
# MPT02           maj, meet (y join z)
# MPT0inf         x meet (y join z)
# MT1k, k >= 2    th_2^k+1, 1
# MT1inf          x join (y meet z), 1
# MPT1k, k >= 3   th_2^k+1 for k ≥ 3,
# MPT12           maj, join (y meet z)
# MPT1inf         x join (y meet z)
# MEET            meet, 0, 1
# MEETP0          meet, 0
# MEETP1          meet, 1
# MEETP           meet
# JOIN            join, 0, 1
# JOINP0          join, 0
# JOINP1          join, 1
# JOINP           join
# D               maj, neg
# DP              maj, x + y + z
# DM              maj
# A               iff, 0
# AD              neg, x + y + z
# AP0             +
# AP1             iff
# AP              x + y + z
# U               neg, 0
# UD              neg
# UM              0, 1
# UP0             0
# UP1             1
# F               { }
#----------------------------------------------------------------------------}}}

# imports {{{1
from random import choice
import ualgebra as UA
#---------------------------------------------------------------------------}}}1

# operations {{{1
def threshold(k, n): # {{{
  """
  Return the threshold function (the UA.Operator version of it) th_k^n defined
  by

    th_k^n(x1, ..., xn) = 1   if #{i | x_i = 1} >= k,
                          0   else
  """
  def thnk(*args, k=k, n=n):
    if args.count(1) >= k:
      return 1
    return 0
  op_thnk = UA.Operation(thnk, n, "th_" + str(k) + "^" + str(n))
  return op_thnk
#----------------------------------------------------------------------------}}}

# 1-ary functions {{{
def const0(x):
  return 0
op_const0 = UA.Operation(const0, 1, "const0")

def const1(x):
  return 1
op_const1 = UA.Operation(const1, 1, "const1")

def neg(x):
  return 1-x
op_neg = UA.Operation(neg, 1, "neg")
#----------------------------------------------------------------------------}}}
# 2-ary functions {{{
def meet(x,y):
  return min(x, y)
op_meet = UA.Operation(meet, 2, "meet")

def join(x,y):
  return max(x, y)
op_join = UA.Operation(join, 2, "join")

def plus(x,y):
  return (x + y) % 2
op_plus = UA.Operation(plus, 2, "plus")

def implies(x,y):
  if x == 1 and y == 0:
    return 0
  return 1
op_implies = UA.Operation(implies, 2, "implies")

def iff(x,y):
  if x == y:
    return 1
  return 0
op_iff = UA.Operation(iff, 2, "iff")

def not_implies(x,y):
  return neg(implies(x,y))
op_not_implies = UA.Operation(not_implies, 2, "not_implies")
#----------------------------------------------------------------------------}}}
# 3-ary functions {{{
op_maj = threshold(2, 3)
op_maj.name = "maj"

def q(x,y,z):   # ternary conditional
  return [z, y][x]
op_q = UA.Operation(q, 3, "q")

def PT0_ternary(x,y,z):
  return meet(x, implies(y,z))
op_PT0_ternary = UA.Operation(PT0_ternary, 3, "PT0_ternary")

def PT1_ternary(x,y,z):
  return join(x, plus(y,z))
op_PT1_ternary = UA.Operation(PT1_ternary, 3, "PT1_ternary")

def MT0_ternary(x,y,z):
  return meet(x, join(y,z))
op_MT0_ternary = UA.Operation(MT0_ternary, 3, "MT0_ternary")

def MT1_ternary(x,y,z):
  return join(x, meet(y,z))
op_MT1_ternary = UA.Operation(MT1_ternary, 3, "MT1_ternary")

def plus_ternary(x,y,z):
  return (x + y + z) % 2
op_3plus = UA.Operation(plus_ternary, 3, "3plus")
#----------------------------------------------------------------------------}}}
#----------------------------------------------------------------------------}}}1

def named_clone(name):  # {{{1
  """
  Returns a basis for the named clone.
  Subscripts and superscripts are given inline, e.g. T1^10 = "T110".
  Infinity is indicated with "inf", e.g. "T0inf"
  """

  if name == "T":
    return [ op_join, op_neg ]
  if name == "P0":
    return [ op_join, op_plus ]
  if name == "P1":
    return [ op_meet, op_implies ]
  if name == "P":
    return [ op_q ]
  if name == "M":
    return [ op_const0, op_const1, op_meet, op_join ]
  if name == "MP0":
    return [ op_const0, op_meet, op_join ]
  if name == "MP1":
    return [ op_const1, op_meet, op_join ]
  if name == "MP":
    return [ op_meet, op_join ]
  if name == "MEET":
    return [ op_const0, op_const1, op_meet ]
  if name == "MEETP0":
    return [ op_const0, op_meet ]
  if name == "MEETP1":
    return [ op_const1, op_meet ]
  if name == "MEETP":
    return [ op_meet ]
  if name == "JOIN":
    return [ op_const0, op_const1, op_join ]
  if name == "JOINP0":
    return [ op_const0, op_join ]
  if name == "JOINP1":
    return [ op_const1, op_join ]
  if name == "JOINP":
    return [ op_join ]
  if name == "D":
    return [ op_neg, op_maj ]
  if name == "DP":
    return [ op_maj, op_3plus ]
  if name == "DM":
    return [ op_maj ]
  if name == "A":
    return [ op_const0, op_iff ]
  if name == "AD":
    return [ op_neg, op_3plus ]
  if name == "AP0":
    return [ op_plus ]
  if name == "AP1":
    return [ op_iff ]
  if name == "AP":
    return [ op_3plus ]
  if name == "U":
    return [ op_const0, op_neg ]
  if name == "UD":
    return [ op_neg ]
  if name == "UM":
    return [ op_const0, op_const1 ]
  if name == "UP0":
    return [ op_const0 ]
  if name == "UP1":
    return [ op_const1 ]
  if name == "F":
    return [ ]

  # infinite families {{{
  if name[:2] == "T0":
    R = [ op_not_implies ]
    if name[2:] != "inf":
      k = int(name[2:])
      if k < 2: raise RuntimeError("k must be >= 2.")
      R.append( threshold(k, k+1) )
    return R

  if name[:3] == "PT0":
    R = [ op_PT0_ternary ]
    if name[3:] != "inf":
      k = int(name[3:])
      if k < 2: raise RuntimeError("k must be >= 2.")
      R.append( threshold(k, k+1) )
    return R

  if name[:2] == "T1":
    R = [ op_implies ]
    if name[2:] != "inf":
      k = int(name[2:])
      if k < 2: raise RuntimeError("k must be >= 2.")
      R.append( threshold(2, k+1) )
    return R

  if name[:3] == "PT1":
    R = [ op_PT1_ternary ]
    if name[3:] != "inf":
      k = int(name[3:])
      if k < 2: raise RuntimeError("k must be >= 2.")
      R.append( threshold(2, k+1) )
    return R

  if name[:3] == "MT0":
    R = [ op_const0 ]
    if name[3:] != "inf":
      k = int(name[3:])
      if k < 2: raise RuntimeError("k must be >= 2.")
      R.append( threshold(k, k+1) )
    else:
      R.append( op_MT0_ternary )
    return R

  if name[:4] == "MPT0":
    if name[4:] != "inf":
      k = int(name[4:])
      if k < 2: raise RuntimeError("k must be >= 2.")
      if k == 2:
        return [ op_maj, op_MT0_ternary ]
      else:
        return [ threshold(k, k+1) ]
    else:
      return [ op_MT0_ternary ]

  if name[:3] == "MT1":
    R = [ op_const1 ]
    if name[3:] != "inf":
      k = int(name[3:])
      if k < 2: raise RuntimeError("k must be >= 2.")
      R.append( threshold(2, k+1) )
    else:
      R.append( op_MT1_ternary )
    return R

  if name[:4] == "MPT1":
    if name[4:] != "inf":
      k = int(name[4:])
      if k < 2: raise RuntimeError("k must be >= 2.")
      if k == 2:
        return [ op_maj, op_MT1_ternary ]
      else:
        return [ threshold(2, k+1) ]
    else:
      return [ op_MT1_ternary ]
  #--------------------------------------------------------------------------}}}
#----------------------------------------------------------------------------}}}1

def rand_clone(max_k): # {{{
  Names = [
    "T"    , "P0"    , "P1"    , "P"    , "M"   , "MP0"   , "MP1"   , "MP"   ,
    "MEET" , "MEETP0", "MEETP1", "MEETP", "JOIN", "JOINP0", "JOINP1", "JOINP",
    "D"    , "DP"    , "DM"   , "A"   , "AD"    , "AP0"   , "AP1"   , "AP"   ,
    "U"    , "UD"    , "UM"   , "UP0" , "UP1"   , "F"
    ] + [ "T0" + str(k)   for k in range(2,max_k+1) ] + [ "T0inf" ]   + \
        [ "PT0" + str(k)  for k in range(2,max_k+1) ] + [ "PT0inf" ]  + \
        [ "T1" + str(k)   for k in range(2,max_k+1) ] + [ "T1inf" ]   + \
        [ "PT1" + str(k)  for k in range(2,max_k+1) ] + [ "PT1inf" ]  + \
        [ "MT0" + str(k)  for k in range(2,max_k+1) ] + [ "MT0inf" ]  + \
        [ "MPT0" + str(k) for k in range(2,max_k+1) ] + [ "MPT0inf" ] + \
        [ "MT1" + str(k)  for k in range(2,max_k+1) ] + [ "MT1inf" ]  + \
        [ "MPT1" + str(k) for k in range(2,max_k+1) ] + [ "MPT1inf" ]
  name = choice(Names)
  return name, named_clone(name)
#----------------------------------------------------------------------------}}}

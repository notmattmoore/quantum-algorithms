# Useful printing algorithms
# Author: Taylor Walenczyk
# Last updated: 08.30.2019

# Prints dictionaries nicely
# In:   d, the dictionary to print; level, the level in the printing process
def pprint_dict(d, level=0):
    ret = with_tabs(level,'{\n')
    for key in d:
        if type(d[key]) is dict:
            ret += with_tabs(level+1,key+':')
            ret += pprint_dict(d[key],level+1)
        else:
            ret += with_tabs(level+1,key+':'+d[key]+',\n')
    ret += with_tabs(level,'}\n')
    return ret

# Prints tabs
# In:   n, the number of tabs to print
def with_tabs(n, inp):
    return '\t'*n+inp


# R. Jay Mashl <rmashl@wustl.edu>

import re

# Harmonize, given gene g and variant v from source src
def harmonize_maf( myvar ):

    muts = myvar.split(';')   # in case of MNPs
    outlist = []
    for m in muts:
        outlist.append( m.strip().replace('X','*').replace('p.*','p.X').replace('p.','') )
    return str(';'.join( outlist ))

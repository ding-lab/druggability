
# R. Jay Mashl <rmashl@wustl.edu>

import re

# Harmonize, given gene g and variant v from source src
def harmonize_maf( myvar ):

    muts = myvar.split(';')   # in case of MNPs
    outlist = []
    for m in muts:
        tmp = m.strip().replace('X','*').replace('p.*','p.X').replace('p.','').replace('%3D','=')
        if re.search('=$',tmp):
            tmp = fix_synonymous( tmp )
        outlist.append( tmp )
    return str(';'.join( outlist ))

def fix_synonymous( myvar ):
    m = re.search(r'([A-Y])(\d+)=$', myvar)
    if m is not None:
        return ''.join([ m[1], m[2], m[1] ])
    else:
        print('# ERROR: unexpected synonymous variant format')

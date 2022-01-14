
# R. Jay Mashl <rmashl@wustl.edu>

import re

# Harmonize, given gene g and variant v from source src
def harmonize_variant( g, v, src ):

    if src == 'civic':
        if v=='BCR-ABL' or re.search(r'BCR-ABL\ ', v):
            v = 'BCR-ABL1'

    return v

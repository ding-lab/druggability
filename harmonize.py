
# R. Jay Mashl <rmashl@wustl.edu>

import re
import config
from enums import *

DEBUG_2=config.DEBUG_2

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


def harmonize_maf_2( myvar, gene, Fasta ):
    '''
    Additional round of reformatting due to non-HGVS convention.
    In particular, position ranges do not include the single-letter amino acids.
    We can fix 'delins' variants, which come denoted  as '>', and some insertions we've seen include *{aa} and **.
    '''
    if re.search(r'>', myvar):
        a = myvar.split('>')
        if len(a) != 2:
            print('ERROR: unexpected split of delins variant %s in union maf' % ( myvar ))
            return myvar
        ins_aa = a[1]

        b = a[0].split('_')
        if len(b) != 2:
            print('ERROR: no range reported in delins variant %s in union maf' % ( myvar ))
            return myvar
        start_pos  = int(b[0])
        c          = re.search(r'^(\d+)([A-Y]{2,})$', b[1])
        end_pos    = int(c[1])
        ref_aa_str = c[2]
        start_aa   = ref_aa_str[ 0]
        end_aa     = ref_aa_str[-1]

        new_myvar = start_aa + str(start_pos) + '_' + end_aa + str(end_pos) + 'delins' + ins_aa
        if DEBUG_2:
            print('# union maf: reformat %s -> %s' % ( myvar, new_myvar ))
        return new_myvar

    if re.search(r'del$', myvar):
        a = myvar[ 0 : len(myvar) - 3 ]
        b = re.search(r'([A-Y]+)([0-9]+)', a)
        if b is None:
            print('ERROR: unable to parse del variant %s in union maf' % (myvar))
            return myvar
        if len(b[1]) == 1:
            new_myvar = myvar
        else:
            start_aa  = b[1][0]
            start_pos = int(b[2])
            end_aa    = b[1][-1]
            end_pos   = start_pos + len(b[1]) - 1
            new_myvar = start_aa + str(start_pos) + '_' + end_aa + str(end_pos) + 'del'

        if DEBUG_2:
            print('# union maf: reformat %s -> %s' % ( myvar, new_myvar ))
        return new_myvar

    if re.search(r'ins', myvar):
        b = re.search(r'^(\d+)(_)(\d+)(ins)(\D*)$', myvar)
        if b is None:
            print('ERROR: unable to parse ins variant %s in union maf' % (myvar))
            return myvar

        if re.search(r'ins$', myvar):   # check to provide blank due to optional
            b.append('')
        start_pos = int(b[1])
        end_pos   = int(b[3])
        ins_seq   = b[5]

        # Use amino acid lookup to reformat
        if gene not in Fasta.keys():
            print('ERROR: unable to reformat variant %s in union maf: gene %s does not exist among fasta sequences' % (myvar, gene))
            return myvar
        start_aa  = Fasta[gene]['fasta'][ start_pos - 1 ]
        end_aa    = Fasta[gene]['fasta'][ end_pos   - 1 ]
        new_myvar = start_aa + str(start_pos) + '_' + end_aa + str(end_pos) + 'ins' + ins_seq
        if DEBUG_2:
            print('# union maf: reformat %s -> %s' % ( myvar, new_myvar ))
        return new_myvar

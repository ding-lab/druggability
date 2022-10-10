
# R. Jay Mashl <rmashl@wustl.edu>

import re
import druggability_databases.config as config
import myglobal
from enums import *
import logging

logger = logging.getLogger(__name__)
logger.setLevel(0)

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
        logger.warning('Unexpected synonymous variant format: %s' % ( myvar ))


def harmonize_maf_2( myvar, gene, Fasta ):
    '''
    Additional round of reformatting due to non-HGVS convention.
    In particular, position ranges do not include the single-letter amino acids.
    We can fix 'delins' variants, which come denoted  as '>'
    '''

    if re.search(r'del$', myvar):
        b = re.search(r'(\D+)([0-9]+)del$', myvar)
        if b is None:
            logger.warning('Unable to parse del variant from union maf. Leaving unharmonized as %s' % (myvar))
            return myvar
        if len(b[1]) == 1:
            new_myvar = myvar
        else:
            start_aa  = b[1][0]
            start_pos = int(b[2])
            end_aa    = b[1][-1]
            end_pos   = start_pos + len(b[1]) - 1
            new_myvar = start_aa + str(start_pos) + '_' + end_aa + str(end_pos) + 'del'

        logger.info('Reformatted %s -> %s in union maf' % ( myvar, new_myvar ))
        return new_myvar

    if re.search(r'>', myvar):
        a = myvar.split('>')
        if len(a) != 2:
            logger.warning('Unexpected input format of delins variant from union maf. Leaving unharmonized as %s' % ( myvar ))
            return myvar
        ins_aa = a[1]

        if not re.search( '_', a[0]):
            logger.warning('No range reported in delins variant from union maf. Leaving unharmonized as %s' % ( myvar ))
            return myvar

        b = re.search(r'(\d+)_(\d+)(\D*)', a[0])
        start_pos  = int(b[1])
        end_pos    = int(b[2])
        ref_aa_str = b[3]
        start_aa   = ref_aa_str[ 0]
        end_aa     = ref_aa_str[-1]

        if not len(ref_aa_str):
            logger.warning('Reference allele unspecified in delins variant %s. Leaving unharmonized' % (myvar))
            return myvar
        else:
            new_myvar = start_aa + str(start_pos) + '_' + end_aa + str(end_pos) + 'delins' + ins_aa
            logger.info('Reformatted %s -> %s in union maf' % ( myvar, new_myvar ))
            return new_myvar

    if re.search(r'ins', myvar):
        b = re.search(r'^(\d+)(_)(\d+)(ins)(\D*)$', myvar)
        if b is None:
            logger.warning('Unable to parse ins variant from union maf. Leaving unharmonized as %s' % (myvar))
            return myvar

        if re.search(r'ins$', myvar):   # check to provide blank due to optional
            b.append('')
        start_pos = int(b[1])
        end_pos   = int(b[3])
        ins_seq   = b[5]

        # Use amino acid lookup to reformat
        if gene not in Fasta.keys():
            logger.warning('Unable to reformat variant %s in union maf: gene %s does not exist among fasta sequences' % (myvar, gene))
            return myvar
        try:
            start_aa  = Fasta[gene]['fasta'][ start_pos - 1 ]
        except:
            logger.warning('Fasta start position lookup failed: gene={}, position={}. Gene length is {}. Leaving unharmonized'.format( gene, start_pos - 1, len(Fasta[gene]['fasta']) ))
            return myvar

        try:
            end_aa    = Fasta[gene]['fasta'][ end_pos   - 1 ]
        except:
            logger.warning('Fasta end position lookup failed: gene={}, position={}. Gene length is {}. Leaving unharmonized'.format( gene, end_pos - 1, len(Fasta[gene]['fasta']) ))
            return myvar

        new_myvar = start_aa + str(start_pos) + '_' + end_aa + str(end_pos) + 'ins' + ins_seq
        logger.info('Reformatted %s -> %s in union maf' % ( myvar, new_myvar ))
        return new_myvar

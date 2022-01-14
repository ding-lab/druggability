
# R. Jay Mashl <rmashl@wustl.edu>

import csv, re
import config
from utils import *
from harmonize import *

DEBUG=config.DEBUG

def load_civic(Variants, Genes, VariantAliases):

    tsv_file = open( config.civic_files['variants'])
    read_tsv = csv.reader(tsv_file, delimiter='\t')
    bReadHeader = True
    hdr_civic_variants = []    # field names

    num_vars_read=0
    num_vars_ignored=0

    for row in read_tsv:
        fields = [ s.strip() for s in row ]
        if bReadHeader:
            hdr_civic_variants = list(fields)
            bReadHeader = False
            continue

        num_vars_read += 1
        variant_id = 'civic:' + fields[0]
        gene = fields[2]
        my_variant = fields[4]
        ignore_status = False

        # Determine whether to ignore the alteration
        if my_variant.lower() in ['phosphorylation', 'mutation', 'underexpression', 'overexpression', 'amplification', 'copy number variation']   or \
           re.search(r'transcript_ablation', fields[19]) or \
           re.search(r'exon', my_variant.lower(), re.IGNORECASE):
            ignore_status = True
        if ignore_status == True:
            num_vars_ignored += 1
            if DEBUG:
                print(' '.join( ['# IGNORED:'] + list( fields[x] for x in (2,4,19) ) ))
            continue

        # Harmonize and categorize variant
        my_variant = harmonize_variant( gene, my_variant, 'civic')
        our_variant_category = 'undeclared'
        tmp_variant_types_list = fields[19].split(',')
        if len(intersection( tmp_variant_types_list,  ['transcript_fusion','gene_fusion'])) > 0:
            our_variant_category = 'fusion'
        else:
            our_variant_category = 'mutation'


        # Accepted variant
        if DEBUG:
            print(' '.join( ['# Accepted:', variant_id, gene, my_variant, our_variant_category]))

            # Store variant
        Variants[variant_id] = {
            'gene':               gene,
            'variant':            my_variant,
            'variant_types_list': tmp_variant_types_list,
            'chrom' :             fields[ 7],
            'pos0' :              fields[ 8],
            'pos1' :              fields[ 9],
            'ensemble_version' :  fields[13],
            'ref_build' :         fields[14],
            'variant_aliases':    fields[25],
            'evidence_list':      [],
            'our_variant_category': our_variant_category,
        }

        # Add variant to list for its gene
        if gene not in Genes.keys():
            Genes[gene] = []

        Genes[gene].append( variant_id )

        # Add variant aliases
        for a in Variants[variant_id]['variant_aliases'].split(','):
            a = a.strip()
            #  Check for rs entry
            if re.search( r'^rs[1-9]+', a, re.IGNORECASE):
                a = a.lower()
                if a not in VariantAliases:
                    VariantAliases[a] = []
                if variant_id not in VariantAliases[a]:
                    VariantAliases[a].append( variant_id )

    tsv_file.close()
    if DEBUG:
        print('# num civic variant records read: ' + str(num_vars_read))
        print('# num civic variant records ignored: ' + str(num_vars_ignored))
    print( Variants.keys() )



def load_oncokb(Variants, Genes, VariantAliases):

    tsv_file = open( config.oncokb_files['variants'])
    read_tsv = csv.reader(tsv_file, delimiter='\t')
    bReadHeader = True
    hdr_oncokb_variants = []    # field names

    num_vars_read=0

    for row in read_tsv:
        fields = [ s.strip() for s in row ]
        if fields[0]=="Gene" and bReadHeader:
            hdr_oncokb_variants = list(fields)
            bReadHeader = False
            continue

        num_vars_read += 1
        gene = fields[0]
        variant = fields[1]

        variant_id = 'oncokb:' + str(num_vars_read)

        # Our labeling of variants
        our_variant_category = 'undeclared'
        if re.search( r'fusion', variant, re.IGNORECASE):
            our_variant_category = 'fusion'
        else:
            our_variant_category = 'mutation'

        Variants[variant_id] = {
            'gene':               gene,
            'variant':            variant,
            'variant_types_list': '',
            'chrom' :             '',
            'pos0' :              '',
            'pos1' :              '',
            'ensemble_version' :  '',
            'ref_build' :         '',
            'variant_aliases':    '',
            'evidence_list':      [],
            'our_variant_category': our_variant_category,
        }

        # Add variant to list for its gene
        if gene not in Genes.keys():
            Genes[gene] = []

        Genes[gene].append( variant_id )

    tsv_file.close()
    if DEBUG:
        print('# num oncokb variant records read: ' + str(num_vars_read))


def load_civic_evidence( Evidence, Variants ):
    tsv_file = open( config.civic_files['evidence'])
    read_tsv = csv.reader(tsv_file, delimiter='\t')
    hdr_civic_evidence = []    # field names

    for row in read_tsv:
        fields = [ s.strip() for s in row ]

        # Header
        if fields[0] == "gene":
            hdr_civic_evidence = list(fields)
            continue

        evidence_id = 'civic:' + fields[20]
        variant_id  = 'civic:' + fields[21]

        Evidence[evidence_id] = {
            'source_db' :  'civic',

            'oncogenicity':     '-',
            'mutation_effect':  '-',

            'disease' :                    fields[ 3],
            'drugs_list_string' :          fields[ 6] if len(fields[6]) else 'No Data',
            'evidence_type' :              fields[ 8],
            'evidence_direction' :         fields[ 9],
            'evidence_level':              fields[10],
            'clinical_significance':       fields[11],
            'clinical_trials_list_string': fields[17],
            'citations':                   [],
            }

        # Account for possible multiple sources in their record
        if fields[14]:
            Evidence[evidence_id]['citations'].append( {'source': fields[14], 'citation_id': fields[13]} )
        if fields[15]:
            Evidence[evidence_id]['citations'].append( {'source': 'ASCO',     'citation_id': fields[15]} )

        # Crosslink tables
        if variant_id in Variants:
            Variants[variant_id]['evidence_list'].append( evidence_id )

    tsv_file.close()


def load_oncokb_evidence( Evidence, Variants ):
    # Reopen alterations file to load more data
    tsv_file = open( config.oncokb_files['variants'])
    read_tsv = csv.reader(tsv_file, delimiter='\t')
    hdr_oncokb_evidence = []    # field names

    num_vars_read=0
    num_genes_found=0
    for row in read_tsv:
        fields = [ s.strip() for s in row ]
        # Header
        if fields[0] == "Gene":
            hdr_oncokb_evidence = list(fields)
            continue

        num_vars_read += 1
        variant_id = 'oncokb:' + str(num_vars_read)
        evidence_id = 'oncokb:' + str( num_vars_read)

        Evidence[evidence_id] = {
            'source_db':  'oncokb',

            'oncogenicity':     fields[2],
            'mutation_effect':  fields[3],

            'disease' :    '-',
            'drugs_list_string' : '-',
            'evidence_type' :     '-',
            'evidence_direction' : '-',
            'evidence_level':      '-',
            'clinical_significance':  '-',
            'clinical_trials_list_string': '-',
            'citations':                   [],
            }

        # Crosslink tables if alteration was accepted originally
        if variant_id in Variants:
            Variants[variant_id]['evidence_list'].append( evidence_id )

    tsv_file.close()

def load_oncokb_therapeutics(Evidence, Variants, Genes):
    tsv_file = open( config.oncokb_files['therapeutics'])
    read_tsv = csv.reader(tsv_file, delimiter='\t')

    for row in read_tsv:
        fields = [ s.strip() for s in row ]
        # Header
        if fields[0] == "Gene":
            continue

        gene = fields[0]
        level = fields[1]
        alteration = fields[2]
        cancer_type = fields[3]
        drug = fields[4]

        if gene not in Genes.keys():
            continue
        # Find the right variant
        for v_id in Genes[gene]:
            if Variants[v_id]['variant'] == alteration:
                for e_id in Variants[v_id]['evidence_list']:
                    if Evidence[e_id]['source_db'] == 'oncokb':
                        Evidence[e_id]['disease']           = cancer_type
                        Evidence[e_id]['drugs_list_string'] = drug
                        Evidence[e_id]['evidence_level']   = level

    tsv_file.close()

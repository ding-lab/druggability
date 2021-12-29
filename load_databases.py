
# R. Jay Mashl <rmashl@wustl.edu>

import csv, re
import config
from utils import *

DEBUG=config.DEBUG

def load_civic(Variants, Genes, VariantAliases):

    tsv_file = open( config.civic_files['variants'])
    read_tsv = csv.reader(tsv_file, delimiter='\t')
    bReadHeader = True
    hdr_civic_variants = []    # field names

    num_vars_read=0

    for row in read_tsv:
        fields = [ s.strip() for s in row ]
        if bReadHeader:
            hdr_civic_variants = list(fields)
            bReadHeader = False
            continue

        num_vars_read += 1
        gene = fields[2]

        variant_id = 'civic:' + fields[0]

        # Screen input and categorize variant
        tmp_variant_types_list = fields[19].split(',')

        # Current: keep simple variants not involved in fusions
        if 'missense_variant' in tmp_variant_types_list:
            if len(intersection( tmp_variant_types_list,  ['transcript_fusion'])) > 0:
                continue

        Variants[variant_id] = {
            'gene':               gene,
            'variant':            fields[ 4],
            'variant_types_list': tmp_variant_types_list,
            'chrom' :             fields[ 7],
            'pos0' :              fields[ 8],
            'pos1' :              fields[ 9],
            'ensemble_version' :  fields[13],
            'ref_build' :         fields[14],
            'variant_aliases':    fields[25],
            'evidence_list':      [],
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

        variant_id = 'oncokb:' + str(num_vars_read)
        Variants[variant_id] = {
            'gene':               gene,
            'variant':            fields[ 1],
            'variant_types_list': '',
            'chrom' :             '',
            'pos0' :              '',
            'pos1' :              '',
            'ensemble_version' :  '',
            'ref_build' :         '',
            'variant_aliases':    '',

            'evidence_list':      [],
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

        # Crosslink tables
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

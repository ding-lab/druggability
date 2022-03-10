
# R. Jay Mashl <rmashl@wustl.edu>

import csv, re
import config
from utils import *
from enums import *
import gzip
import logging

DEBUG=config.DEBUG
DEBUG_2=config.DEBUG_2


def load_civic(Variants, Genes, VariantAliases):
    # Load preprocessed variants
    tsv_file = open( config.civic_files['variants_preprocessed'])
    read_tsv = csv.reader(tsv_file, delimiter='\t')
    bReadHeader = True
    num_vars_read=0
    num_vars_ignored=0

    for row in read_tsv:
        fields = [ s.strip() for s in row ]
        if len(fields) != 11:
            abort_run('check column count in preprocessed variants file')
        if bReadHeader:
            if fields[0] == 'variant_id':
                bReadHeader = False
            continue
        num_vars_read += 1
        variant_id         = 'civic:{myid}'.format( myid=fields[0] )
        variant            = fields[3]
        variant_types_str  = fields[4]
        num_conditions     = int(fields[6])   # number of conditions to consider

        # check to ignore
        if num_conditions == 0:    # ignore entry
            num_vars_ignored += 1
            if DEBUG_2:
                logging.info(' '.join( ['preprocessed IGNORED:', variant_id, fields[1],variant,variant_types_str] ))
            continue

        if re.search(r'amplification|copy number|exon|expression|loss|microsatellite|mutation|nuclear|phosphorylation|promoter|splice|tandem|truncating|wildtype|domain', variant.lower(), re.IGNORECASE):
            num_vars_ignored += 1
            if DEBUG_2:
                logging.info(' '.join( ['preprocessed IGNORED(variant name):', variant_id, fields[1],variant,variant_types_str] ))
            continue

        # Classify variant
        main_variant_class = UNDECLARED
        fusion_gene_set = []
        x0, x1 = [0, 0]

        vartypes_list = [ s.strip() for s in variant_types_str.split(',') ]

        if len(intersection( vartypes_list,  ['transcript_fusion','gene_fusion'])) > 0:
            main_variant_class = FUSION

            # Reformat to gene1--gene2
            variant = variant.replace('-','--')

            # Handle exceptions

            # Extract genes
            fusion_gene_set = variant.split('--')

        else:
            main_variant_class = MUTATION

            # Find position range (integral values) for AA change
            x0, x1 = get_pos_range( variant )

        if DEBUG_2:
            logging.info( '(civic) ' + variant + ' is range ' + str(x0) + ' - ' + str(x1))

        # Instantiate record
        Variants[variant_id] = {
            'gene'                : fields[ 1],
            'variant'             : variant,        # can contain manually standardized aachange
            'vartypes_list'       : vartypes_list,

            'pp_comment'          : fields[ 5],
            'pp_conditions'       : num_conditions,
            'pp_condition2'       : fields[ 7],
            'pp_condition2_value' : fields[ 8],
            'pp_condition3'       : fields[ 9],
            'pp_condition3_value' : '' if fields[10] == '-' else fields[10],
            'fusion_gene_set'     : fusion_gene_set,
            'main_variant_class'  : main_variant_class,

            'chrom'               : '',
            'pos0'                : '',
            'pos1'                : '',
            'ref'                 : '',
            'alt'                 : '',
            'ensemble_version'    : '',
            'ref_build'           : '',
            'variant_aliases'     : '',
            'evidence_list'       : [],
            'gdnachange'          : '',
            'gdnacoords'          : '',

            'prot_ref_start_pos'  : x0,  # start, end positions in protein target by maf mutation
            'prot_ref_end_pos'    : x1,

            'comment_liftover'    : '',
            'chrom_liftover'      : '',
            'start_liftover'      : '',
            'stop_liftover'       : '',
            'ref_build_liftover'  : '',
            'gdnachange_liftover' : '',
            'gdnacoords_liftover' : '',
        }

    tsv_file.close()
    if DEBUG:
        logging.info('num preprocessed civic variant records read: ' + str(num_vars_read))
        logging.info('num preprocessed civic variant records ignored: ' + str(num_vars_ignored))


    # Load (modified) native file to get more information
    tsv_file = open( config.civic_files['variants'])
    read_tsv = csv.reader(tsv_file, delimiter='\t')
    bReadHeader = True
    hdr_civic_variants = []    # field names

    num_vars_read    = 0
    num_vars_ignored = 0

    for row in read_tsv:
        fields = [ s.strip() for s in row ]
        if bReadHeader:
            hdr_civic_variants = list(fields)
            bReadHeader        = False
            continue

        num_vars_read += 1
        variant_id = 'civic:{myid}'.format( myid = fields[0] )

        # Safeguard against new entries
        if variant_id not in Variants.keys():
            continue

        # Omit liftover failures by looking at chromosome field
        if fields[30] == '-':
            continue

        # Accepted variant
        if DEBUG_2:
            logging.info(' '.join( ['Accepted:', variant_id, fields[1] + ' as ' + Variants[variant_id]['gene'], fields[3]] ))

        # Update information
        chrom = fields[ 7].replace('chr','')
        pos0  = fields[ 8]
        pos1  = fields[ 9]
        Variants[variant_id].update({
            'chrom'              : chrom,
            'pos0'               : pos0,
            'pos1'               : pos1,
            'ref'                : fields[10],
            'alt'                : fields[11],
            'ensemble_version'   : fields[13],
            'ref_build'          : fields[14],
            'variant_aliases'    : fields[25],
            'evidence_list'      : [],
            'comment_liftover'   : fields[29],
            'chrom_liftover'     : fields[30].replace('chr',''),
            'start_liftover'     : fields[31],
            'stop_liftover'      : fields[32],
            'ref_build_liftover' : config.civic_params['ref_build_liftover'],
        })

        # Validate input
        if not len(chrom):
            abort_run('chromosome is missing from data source. Please check the input.')
        else:
            if not( len(pos0) and len(pos1) ):
                abort_run('this data source typically specifies chr,start,stop coordinates. Please check the input.')
            if int(pos1) < int(pos0):
                abort_run('stop coordinate is upstream from start. Please check the input.')

        # Calculate gDNA features
        # ...for convenience, if not overkill
        Variants[variant_id].update({
            'gdnachange': calculate_gdna_change(Variants[variant_id], ''),
            'gdnacoords': calculate_gdna_coords(Variants[variant_id], ''),

            'gdnachange_liftover': calculate_gdna_change(Variants[variant_id], 'use_liftover'),
            'gdnacoords_liftover': calculate_gdna_coords(Variants[variant_id], 'use_liftover'),
        })

        # Add variant to gene list
        gene = Variants[variant_id]['gene']
        if gene not in Genes.keys():
            Genes[gene] = []
        Genes[gene].append( variant_id )

        # Track variant aliases
        for a in Variants[variant_id]['variant_aliases'].split(','):
            a = a.strip()
            #  Check for rs entry
            if re.search( r'^rs[1-9]+', a, re.IGNORECASE):
                a = a.lower()
                if a not in VariantAliases:
                    VariantAliases[a] = []
                list_append( VariantAliases[a], variant_id )

    tsv_file.close()
    if DEBUG:
        logging.info('num civic variant records read: ' + str(num_vars_read))
        logging.info('num civic variant records ignored: ' + str(num_vars_ignored))



def load_oncokb(Variants, Genes, VariantAliases):

    tsv_file = open( config.oncokb_files['variants'])
    read_tsv = csv.reader(tsv_file, delimiter='\t')
    bReadHeader = True
    hdr_oncokb_variants = []    # field names

    num_vars_read=0
    num_vars_ignored=0
    for row in read_tsv:
        fields = [ s.strip() for s in row ]
        if fields[0]=="Gene" and bReadHeader:
            hdr_oncokb_variants = list(fields)
            bReadHeader = False
            continue

        num_vars_read += 1
        gene       = fields[0]
        variant    = fields[1]
        variant_id = 'oncokb:{myid}'.format( myid = str(num_vars_read) )

        # Ignore certain values
        if re.search(r'splice|expression|Dx2|[0-9](mis|mut)$|wildtype|fusions|amplification|truncating|deletion|domain|exon|trunc$|tandem|methylation', variant, re.IGNORECASE):
            num_vars_ignored += 1
            continue

        # Classify variant
        main_variant_class = UNDECLARED
        fusion_gene_set = []
        x0, x1 = [0, 0]

        if re.search( r'fusion', variant, re.IGNORECASE):
            main_variant_class = FUSION

            # Reformat to gene1--gene2
            variant = variant.replace('Fusion', '').strip().replace('-','--')

            # Handle exceptions
            variant = variant.replace('HLA--DRB1','HLA-DRB1')

            # Extract genes
            fusion_gene_set = variant.split('--')

        else:
            main_variant_class = MUTATION

            # Handle exceptions
            if gene in ['EGFR']  and  re.search(r'v[IVX]{1,}', variant):    # exonic isoform
                num_vars_ignored += 1
                continue

            # Find position range (integral values)
            x0, x1 = get_pos_range( variant )

        if DEBUG_2:
            logging.info(' '.join( ['oncokb variant ALLOWED:', variant_id, gene, variant, 'class='+str(main_variant_class)] ))
            logging.info( '(oncokb) ' + variant + ' is range ' + str(x0) + ' - ' + str(x1))


        # Instantiate record
        Variants[variant_id] = {
            'gene'                : gene,
            'variant'             : variant,
            'vartypes_list'       : [],

            'pp_comment'          : '',
            'pp_conditions'       : 1,
            'pp_condition2'       : '',
            'pp_condition2_value' : '',
            'pp_condition3'       : '',
            'pp_condition3_value' : '',
            'fusion_gene_set'     : fusion_gene_set,
            'main_variant_class'  : main_variant_class,

            'chrom'               : '',
            'pos0'                : '',
            'pos1'                : '',
            'ensemble_version'    : '',
            'ref_build'           : '',
            'variant_aliases'     : '',
            'evidence_list'       : [],
            'gdnachange'          : '',
            'gdnacoords'          : '',

            'prot_ref_start_pos'  : x0,  # start, end positions in protein target by maf mutation
            'prot_ref_end_pos'    : x1,

            'comment_liftover'    : '',
            'chrom_liftover'      : '',
            'start_liftover'      : '',
            'stop_liftover'       : '',
            'ref_build_liftover'  : '',
            'gdnachange_liftover' : '',
            'gdnacoords_liftover' : '',
        }

        # Add variant to list for its gene
        if gene not in Genes.keys():
            Genes[gene] = []
        Genes[gene].append( variant_id )

    tsv_file.close()
    if DEBUG:
        logging.info('num oncokb variant records read (initial pass): %s' %  (num_vars_read))
        logging.info('num oncokb variant records ignored (initial pass): %s' % (num_vars_ignored))


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

        evidence_id = 'civic:{myid}'.format( myid = fields[20] )
        variant_id  = 'civic:{myid}'.format( myid = fields[21] )

        Evidence[evidence_id] = {
            'source_db'                  : 'civic',

            'oncogenicity'               : '-',
            'mutation_effect'            : '-',

            'disease'                    : fields[ 3],
            'drugs_list_string'          : fields[ 6] if len(fields[6]) else 'No Data',
            'evidence_type'              : fields[ 8],
            'evidence_direction'         : fields[ 9],
            'evidence_level'             : fields[10],
            'clinical_significance'      : fields[11],
            'clinical_trials_list_string': fields[17],
            'citations'                  : [],
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

    # Reopen alterations file to load evidence
    tsv_file = open( config.oncokb_files['variants'])
    read_tsv = csv.reader(tsv_file, delimiter='\t')
    hdr_oncokb_evidence = []    # field names

    num_vars_read    = 0
    num_vars_ignored = 0
    num_genes_found  = 0

    for row in read_tsv:
        fields = [ s.strip() for s in row ]
        # Header
        if fields[0] == "Gene":
            hdr_oncokb_evidence = list(fields)
            continue

        num_vars_read += 1
        variant_id     = 'oncokb:{myid}'.format( myid = str(num_vars_read) )
        evidence_id    = 'oncokb:{myid}'.format( myid = str(num_vars_read) )

        # proceed only for records retained previously in load_oncokb()
        if  variant_id  not in  Variants:
            num_vars_ignored += 1
            continue

        Evidence[evidence_id] = {
            'source_db'                  : 'oncokb',

            'oncogenicity'               : fields[2],
            'mutation_effect'            : fields[3],

            'disease'                    : '-',
            'drugs_list_string'          : '-',
            'evidence_type'              : '-',
            'evidence_direction'         : '-',
            'evidence_level'             : '-',
            'clinical_significance'      : '-',
            'clinical_trials_list_string': '-',
            'citations'                  : [],
            }

        # Crosslink tables
        Variants[variant_id]['evidence_list'].append( evidence_id )

    tsv_file.close()
    if DEBUG:
        logging.info('num oncokb variant records read (second pass): %s' % (num_vars_read))
        logging.info('num oncokb variant records ignored (second pass): %s' % (num_vars_ignored))


def load_oncokb_therapeutics(Evidence, Variants, Genes):
    tsv_file = open( config.oncokb_files['therapeutics'])
    read_tsv = csv.reader(tsv_file, delimiter='\t')

    for row in read_tsv:
        fields = [ s.strip() for s in row ]
        # Header
        if fields[0] == "Gene":
            continue

        gene        = fields[0]
        level       = fields[1]
        alteration  = fields[2]
        cancer_type = fields[3]
        drug        = fields[4]

        if gene not in Genes.keys():
            continue
        # Find the right variant
        for v_id in Genes[gene]:
            if Variants[v_id]['variant'] == alteration:
                for e_id in Variants[v_id]['evidence_list']:
                    if Evidence[e_id]['source_db'] == 'oncokb':
                        Evidence[e_id]['disease']           = cancer_type
                        Evidence[e_id]['drugs_list_string'] = drug
                        Evidence[e_id]['evidence_level']    = level

    tsv_file.close()


def load_fasta( Fasta ):
    tsv_file = gzip.open( config.uniprot_files['fasta'], mode='rt')
    read_tsv = csv.reader(tsv_file, delimiter='\t')
    bReadHeader = True

    for row in read_tsv:
        fields = [ s.strip() for s in row ]
        if len(fields) != 3:
            abort_run('check column count in preprocessed uniprot fasta file')
        if bReadHeader:
            if fields[0] == 'Gene':
                bReadHeader = False
            continue

        gene, fa_header, fa_seq = fields
        Fasta[ gene ] = {'fasta': fa_seq}
    tsv_file.close()

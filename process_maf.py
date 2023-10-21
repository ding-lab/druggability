
# R. Jay Mashl <rmashl@wustl.edu>

import os, sys, csv, re
import druggability_databases.config as config
import myglobal
from utils import *
from enums import *
from harmonize import *
import logging
from annotate_trials import *

logger = logging.getLogger(__name__)
logger.setLevel(0)

# accommodate large strings
csv.field_size_limit( 131072 * 82 )


# Block to determine input format
def is_maf_header_1( fields ):   # washu
     if ['Chromosome', 'Start_Position', 'End_Position', 'Variant_Type', 'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode', 'HGVSp_Short'] == [ fields[x] for x in [4,5,6,9,15,16,36] ]:
         return True
     else:
         return False

def is_maf_header_2( fields ):   # union
    if ['Chromosome', 'Start_Position', 'End_Position', 'Variant_Type', 'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode', 'Protein_Change'] == [ fields[x] for x in [3,4,5,8,12,13,26] ]:
        return True
    else:
        return False

def is_maf_header_3( fields ):   # basic
    if ['Chromosome', 'Start_Position', 'End_Position', 'Variant_Type', 'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode', ] == [ fields[x] for x in [2,3,4,6,9,10] ] and fields[12] in ['Protein_Change', 'HGVSp_Short', 'HGVSp']:
        return True
    else:
        return False

def get_maf_type( fields ):
    NF = len(fields)
    if NF < 13:
        abort_run('Unrecognized maf format')
    elif NF < 27:
        if is_maf_header_3( fields ):
            return BASIC_MAF
        else:
            abort_run('Unrecognized maf format')
    elif NF < 37:
        if is_maf_header_2( fields ):
            return UNION_MAF
        elif is_maf_header_3( fields ):
            return BASIC_MAF
        else:
            abort_run('Unrecognized maf format')
    else:
        if is_maf_header_1( fields ):
            return WASHU_MAF
        elif is_maf_header_2( fields ):
            return UNION_MAF
        elif is_maf_header_3( fields ):
            return BASIC_MAF
        else:
            abort_run('Unrecognized maf format')


# Process maf (or VCF)
def process_maf( args, Matches, Evidence, Variants, Genes, Fasta, Genes_altered, Trials, Matches_trials, SampleMentioned, call_context, GenesSeenInTrials ):
    #Variant_tracking = dict()   # record variants by sample
    maf_filetype     = UNDECLARED

    bReadHeader      = True
    #bHasSampleMatch   = False
    #SampleMentioned['maf'] = False
    logger.info('In process_maf: Sample mentioned = {}'. format(SampleMentioned['maf']))

    tsv_file = open( args.varInputFile )
    read_tsv = csv.reader(tsv_file, delimiter='\t')

    # Set composite key for tracking matches
    #sample_pair = '{tumor}||{normal}'.format( tumor=args.tumor_name, normal=args.normal_name )
    sample_pair = SAMPLENAME


    # Main
    for row in read_tsv:
        fields = [ s.strip() for s in row ]

        if bReadHeader:
            if row[0] == 'Hugo_Symbol':
                maf_filetype = get_maf_type( fields )
            elif re.search(r'CHROM$', row[0]) and row[1:5] == ['POS','ID','REF','ALT']:
                if row[5] == 'set':
                    maf_filetype = CLINICAL_VCF
                elif row[5:8] == ['QUAL','FILTER','INFO']:
                    maf_filetype = BASIC_VCF

            if maf_filetype != UNDECLARED:
                logger.info( 'maf/vcf format recognized: {}'.format(map_maf_reverse(maf_filetype)) )
                bReadHeader = False

            continue

        if SampleMentioned['maf'] < 0:
            SampleMentioned['maf'] = 0     # There is/are variants listed in input

        if maf_filetype == WASHU_MAF:
            gene         = fields[ 0]
            ref_build    = fields[ 3]
            chrom        = fields[ 4].replace('chr','')
            pos_start    = fields[ 5]
            pos_end      = fields[ 6]
            maf_varclass = fields[ 8]
            vartype      = fields[ 9]   # SNP, INS, DEL, ...

            ref = fields[10]

            # Find an allele that is different: may be a nucleotide(s) or '-'
            for c in [ fields[11], fields[12] ]:
                if c != ref:
                    alt = c

            sample_t     = fields[15]   # here, this is tumor sample
            sample_n     = fields[16]   # here, this is (matched) normal sample
            cdnachange   = fields[34]
            aachange     = fields[36]

        if maf_filetype == UNION_MAF:
            gene         = fields[ 0]
            ref_build    = fields[ 2]
            chrom        = fields[ 3].replace('chr','')
            pos_start    = fields[ 4]
            pos_end      = fields[ 5]
            maf_varclass = fields[ 7]
            vartype      = fields[ 8]   # [SDOTM]NP, INS, DEL

            ref = fields[ 9]

            if ref == fields[10]:
                alt = fields[11] if fields[11] != '-' else ''
            elif ref == fields[11]:
                alt = fields[10] if fields[10] != '-' else ''
            else:
                abort_run('cannot resolve alteration at {gene} {pos}'.format( gene=gene, pos=str(pos_start) ))


            sample_t     = fields[12]   # here, this is tumor sample
            sample_n     = fields[13]   # here, this is (matched) normal sample
            cdnachange   = fields[24]
            aachange     = fields[26]

        if maf_filetype == BASIC_MAF:
            gene         = fields[ 0]
            ref_build    = fields[ 1]
            chrom        = fields[ 2].replace('chr','')
            pos_start    = fields[ 3]
            pos_end      = fields[ 4]
            maf_varclass = fields[ 5]
            vartype      = fields[ 6]   # [SDOTM]NP, INS, DEL
            ref          = fields[ 7]
            alt          = fields[ 8]
            sample_t     = fields[ 9]   # might be same as case id
            sample_n     = fields[10]   # might be same as case id
            cdnachange   = fields[11]
            aachange     = fields[12]

        if maf_filetype == CLINICAL_VCF:
            gene         = fields[15]
            ref_build    = 'ND'   # not in output records
            chrom        = fields[ 0].replace('chr','')
            pos_start    = fields[ 1]
            maf_varclass = fields[14]
            cdnachange   = fields[18].split(':')[1] if re.search(r'^ENST',fields[18]) else ''

            ref = fields[ 3]
            alt = fields[ 4]

            # Calculated values
            RL = len(ref)
            AL = len(alt)
            if RL == AL:   # xNP
                pos_end = str( int(pos_start) + RL - 1 )
                if RL < 4:
                    vartype = [0,'SNP','DNP','TNP'][RL]
                else:
                    vartype = 'MNP'
            elif RL < AL:
                vartype = 'INS'
                pos_end = str( int(pos_start) + 1 )
            else:
                vartype = 'DEL'
                pos_end = str( int(pos_start) + RL - 1 )

            protein_pos = fields[22]
            somatic_aas = fields[23]
            # Use amino acid field as protein position field may be nonempty/not useful
            if len(somatic_aas):
                tmp_aas = somatic_aas.split('/')
                if len(tmp_aas) == 1:
                    aachange = 'p.' + tmp_aas[0] + protein_pos + tmp_aas[0]
                else:
                    aachange = 'p.' + tmp_aas[0] + protein_pos + tmp_aas[1]
            else:
                aachange = ''

            # The pipeline we are accessing does not include samples in its output records so just declare them
            sample_t = args.tumor_name
            sample_n = args.normal_name

        if maf_filetype == BASIC_VCF:
            abort_run('Standard vcf parsing not yet fully implemented')


        # Require samples to match
        if (sample_t != args.tumor_name)  or  (sample_n != args.normal_name):
            continue
        #bHasSampleMatch = True
        SampleMentioned['maf'] = 1

        # Make a list of all genes in the maf
        Genes_altered[ gene ] = 1

        # summarize alteration
        alteration_summary = '\t'.join([ gene, chrom, pos_start, cdnachange, aachange, ref_build, maf_varclass, vartype ])

        # initially look at only vars with AA change in HGVS short format; if blank, it is often a splice site
        if not re.search( r'^p\.', aachange):
            if myglobal.DEBUG_2:
                logger.info( 'maf/vcf record IGNORED: {}' . format(alteration_summary ))
            continue

        # harmonize
        aachange = harmonize_maf( aachange )

        if maf_filetype == UNION_MAF:       # additional treatment
            if re.search(r'>', aachange) or re.search('ins', aachange) or re.search(r'del$', aachange):
                aachange = harmonize_maf_2( aachange, gene, Fasta )

        # As this is a maf, assume no fusions; classification intended for comparison to alteration db
        # (See https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/#protected-maf-file-structure)
        main_variant_class = MUTATION

        # calculate gdna change
        tmp_set = { 'chrom': chrom, 'pos0': pos_start, 'pos1': pos_end, 'ref': ref, 'alt': alt, 'gene': gene, 'main_variant_class': main_variant_class }
        gdnachange = calculate_gdna_change( tmp_set, '' )
        gdnacoords = calculate_gdna_coords( tmp_set, '' )

        # Merge lists of reported rs ids; no such seems available in union_maf
        known_variants = []
        if maf_filetype == WASHU_MAF:
            dbSNP_RS                       = fields[ 13]  #rs#, list of rs#s, or novel
            existing_variation_list_string = fields[ 57]
            maf_variant_id                 = fields[111]
            if dbSNP_RS != '.':
                known_variants.extend( dbSNP_RS.split(','))
            if existing_variation_list_string != '.':
                known_variants.extend( existing_variation_list_string.split(','))
            if maf_variant_id != '.':
                known_variants.extend( maf_variant_id.split(','))
            known_variants = uniquify(  known_variants )

        # Set up storage for tracking matches
        check_alloc_named( Matches, sample_pair, 'match_level' )

        # Check for matches with alterations database
        if gene in Genes.keys():

            # Identify and track matches
            num_full_matches    = 0
            num_partial_matches = 0
            num_unmatched       = 0

            for v_id in Genes[gene]:

                if has_genomic_match( gdnachange, v_id, Variants, 'gdnachange_liftover' ):
                    called_str = 'MUT:{gene} {var}|{genomic_change}|{ref}'.format( gene=gene, var=aachange, genomic_change=gdnachange, ref=ref_build )
                    list_append( Matches[ sample_pair ]['full'], {'v_id': v_id, 'reason': '1.gcoordExact_altExact', 'called': called_str} )
                    continue

                if has_genomic_match( gdnacoords, v_id, Variants, 'gdnacoords_liftover' ):
                    called_str = 'MUT:{gene} {var}|{genomic_change}|{ref}'.format( gene=gene, var=aachange, genomic_change=gdnachange, ref=ref_build )
                    list_append( Matches[ sample_pair ]['partial'], {'v_id': v_id, 'reason': '2.gcoordExact_altDiff', 'called': called_str} )
                    continue

                if has_genomic_overlap( tmp_set, v_id, Variants ):
                    called_str = 'MUT:{gene} {var}|{genomic_change}|{ref}'.format( gene=gene, var=aachange, genomic_change=gdnachange, ref=ref_build )
                    list_append( Matches[ sample_pair ]['partial'], {'v_id': v_id, 'reason': '3.gcoordOverlap', 'called': called_str} )
                    continue

                if Variants[v_id]['main_variant_class'] == MUTATION:
                    if Variants[v_id]['prot_ref_start_pos'] > 0:
                        if get_aachange_overlap_length( aachange, v_id, Variants ) > 0:
                            called_str = 'MUT:{gene} {var}|{genomic_change}|{ref}'.format( gene=gene, var=aachange, genomic_change=gdnachange, ref=ref_build )
                            list_append( Matches[ sample_pair ]['partial'], {'v_id': v_id, 'reason': '4.pcoordOverlap', 'called': called_str} )
                            continue

        # Check for potentially relevant trials
        if args.annotate_trials:

            if re.search(r'[SDOTM]NP', vartype):
                if gene in Trials.keys():
                    if call_context in Trials[ gene ][ MUTATION ].keys():
                        if len( Trials[ gene ][ MUTATION ][ call_context ].keys() ):
                            GenesSeenInTrials[ MUTATION ].append( gene )

                            # check for nonspecific match
                            if 'any' in Trials[ gene ][ MUTATION ][ call_context ].keys():
                                check_alloc_named( Matches_trials[ sample_pair ][ MUTATION ], gene, 'list' )
                                Matches_trials[ sample_pair ][ MUTATION ][ gene ].extend( Trials[ gene ][ MUTATION ][ call_context ][ 'any' ] )
                                if myglobal.DEBUG:
                                    logger.info('{gene} matched -any- trials for xNP' . format(gene=gene))

                            # check for aachanges "as is" or with 'any' alternative allele tag
                            p = aachange
                            m = re.search(r'([A-Z]\d+)(.*)', p )
                            any_probe = m.groups()[0] + ':any:'

                            if myglobal.DEBUG:
                                logger.info('we will check for vartype={x1}, gene={x2}, pos={x3}, any={x4}' . format(x1=vartype, x2=gene, x3=p, x4=any_probe))

                            if  p         in Trials[ gene ][ MUTATION ][ call_context ].keys():
                                check_alloc_named( Matches_trials[ sample_pair ][ MUTATION ], gene, 'list' )
                                Matches_trials[ sample_pair ][ MUTATION ][ gene ].extend( Trials[ gene ][ MUTATION ][ call_context ][ p ] )
                                if myglobal.DEBUG:
                                    logger.info('...matched aachange as {}' . format(p))

                            if any_probe in Trials[ gene ][ MUTATION ][ call_context ].keys():
                                check_alloc_named( Matches_trials[ sample_pair ][ MUTATION ], gene, 'list' )
                                Matches_trials[ sample_pair ][ MUTATION ][ gene ].extend( Trials[ gene ][ MUTATION ][ call_context ][ any_probe ] )
                                if myglobal.DEBUG:
                                    logger.info('...matched aachange as any probe {}' . format(any_probe))

                            # check for genomic coordinate
                            gc_probe = ':'.join([ chrom, pos_start, pos_end, ref, alt ])
                            if gc_probe in Trials[ gene ][ MUTATION ][ call_context ].keys():
                                check_alloc_named( Matches_trials[ sample_pair ][ MUTATION ], gene, 'list' )
                                Matches_trials[ sample_pair ][ MUTATION ][ gene ].extend( Trials[ gene ][ MUTATION ][ call_context ][ gc_probe ] )
                                if myglobal.DEBUG:
                                    logger.info('...matched aachange as any probe {}' . format(gc_probe))


            elif vartype in ['INS', 'DEL']:
                vt = INSERTION if vartype=='INS' else DELETION
                if gene in Trials.keys():
                    if call_context in Trials[ gene ][ vt ].keys():
                        if len( Trials[ gene ][ vt ][ call_context ].keys() ):
                            GenesSeenInTrials[ vt ].append( gene )

                            # check for nonspecific match
                            if 'any' in Trials[ gene ][ vt ][ call_context ].keys():
                                check_alloc_named( Matches_trials[ sample_pair ][ vt ], gene, 'list' )
                                Matches_trials[ sample_pair ][ vt ][ gene ].extend( Trials[ gene ][ vt ][ call_context ][ 'any' ] )
                                if myglobal.DEBUG:
                                    logger.info('{gene} matched -any- trials for {vartype}' . format(gene=gene, vartype=map_mut_reverse(vt)))

                            # check for aachanges "as is"
                            if vt == INSERTION:
                                if alt in Trials[ gene ][ vt ][ call_context ].keys():
                                    check_alloc_named( Matches_trials[ sample_pair ][ vt ], gene, 'list' )
                                    Matches_trials[ sample_pair ][ vt ][ gene ].extend( Trials[ gene ][ vt ][ call_context ][ alt ] )
                                    if myglobal.DEBUG:
                                        logger.info('{gene} matched trials for {vartype} with alteration {alteration}' . format(gene=gene, vartype=map_mut_reverse(vt), alteration=alt))

                            elif vt == DELETION:
                                if ref in Trials[ gene ][ vt ][ call_context ].keys():
                                    check_alloc_named( Matches_trials[ sample_pair ][ vt ], gene, 'list' )
                                    Matches_trials[ sample_pair ][ vt ][ gene ].extend( Trials[ gene ][ vt ][ call_context ][ ref ] )
                                    if myglobal.DEBUG:
                                        logger.info('{gene} matched trials for {vartype} with alteration {alteration}' . format(gene=gene, vartype=map_mut_reverse(vt), alteration=ref))

                            else:
                                pass

            else:
                logger.warning('Unexpected vartype {} in input variant file' . format(vartype) )

    tsv_file.close()

    if bReadHeader:
        abort_run('Variant file {} is missing the header line' . format( args.varInputFile ))

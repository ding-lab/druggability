
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


# Process somatic mafs
def process_maf( args, Matches, Evidence, Variants, Genes, Fasta, Genes_altered, Trials, Matches_trials, call_mode ):
    inputFile        = args.variant_file
    Variant_tracking = dict()   # record variants by sample
    maf_filetype     = UNDECLARED

    bReadHeader      = True
    bHasSampleMatch   = False

    tsv_file = open( inputFile )
    read_tsv = csv.reader(tsv_file, delimiter='\t')

    # Set composite key for tracking matches
    sample_pair = '{tumor}||{normal}'.format( tumor=args.tumor_name, normal=args.normal_name )

    # Track genes lists having matches to trials
    GenesSeenInTrials = dict()
    for alt_type in [ MUTATION, INSERTION, DELETION ]:
        GenesSeenInTrials[ alt_type ] = []

    # Set up storage
    if len(args.annotate_trials):
        check_alloc_named( Matches_trials, sample_pair, 'dict')
        for vc in VARIANT_CLASSES:
            check_alloc_named( Matches_trials[ sample_pair ], vc, 'dict' )

    # Main
    for row in read_tsv:
        fields = [ s.strip() for s in row ]

        if bReadHeader:
            if row[0] == 'Hugo_Symbol':
                maf_filetype = get_maf_type( fields )
                logger.info( 'maf format recognized: {}'.format(map_maf_reverse(maf_filetype)) )
                bReadHeader = False
            continue

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


        # Require samples to match
        if (sample_t != args.tumor_name)  or  (sample_n != args.normal_name):
            continue
        bHasSampleMatch = True

        # Track altered genes seen
        Genes_altered[ gene ] = 1

        # summarize alteration
        alteration_summary = '\t'.join([ gene, chrom, pos_start, cdnachange, aachange, ref_build, maf_varclass, vartype ])

        # initially look at only vars with AA change in HGVS short format; if blank, it is often a splice site
        if not re.search( r'^p\.', aachange):
            if myglobal.DEBUG_2:
                logger.info( 'maf record IGNORED: {}' . format(alteration_summary ))
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

        if sample_pair not in Variant_tracking.keys():
            Variant_tracking[ sample_pair ] = dict()
        Variant_tracking[ sample_pair ][alteration_summary] = dict( total_evidence_count=0, v_id_list=[] )

        # Check for matches with alterations database
        if gene in Genes.keys():

            # Identify and track matches
            num_full_matches    = 0
            num_partial_matches = 0
            num_unmatched       = 0

            for v_id in Genes[gene]:

                if has_genomic_match( gdnachange, v_id, Variants, 'gdnachange_liftover' ):
                    called_str = '{gene} {var}|{genomic_change}|{ref}'.format( gene=gene, var=aachange, genomic_change=gdnachange, ref=ref_build )
                    list_append( Matches[ sample_pair ]['full'], {'v_id': v_id, 'reason': '1.gcoordExact_altExact', 'called': called_str} )
                    Variant_tracking[sample_pair][alteration_summary]['v_id_list'].append( v_id )
                    Variant_tracking[sample_pair][alteration_summary]['total_evidence_count'] += len(Variants[v_id]['evidence_list'])
                    continue

                if has_genomic_match( gdnacoords, v_id, Variants, 'gdnacoords_liftover' ):
                    called_str = '{gene} {var}|{genomic_change}|{ref}'.format( gene=gene, var=aachange, genomic_change=gdnachange, ref=ref_build )
                    list_append( Matches[ sample_pair ]['partial'], {'v_id': v_id, 'reason': '2.gcoordExact_altDiff', 'called': called_str} )
                    Variant_tracking[sample_pair][alteration_summary]['v_id_list'].append( v_id )
                    Variant_tracking[sample_pair][alteration_summary]['total_evidence_count'] += len(Variants[v_id]['evidence_list'])
                    continue

                if has_genomic_overlap( tmp_set, v_id, Variants ):
                    called_str = '{gene} {var}|{genomic_change}|{ref}'.format( gene=gene, var=aachange, genomic_change=gdnachange, ref=ref_build )
                    list_append( Matches[ sample_pair ]['partial'], {'v_id': v_id, 'reason': '3.gcoordOverlap', 'called': called_str} )
                    Variant_tracking[sample_pair][alteration_summary]['v_id_list'].append( v_id )
                    Variant_tracking[sample_pair][alteration_summary]['total_evidence_count'] += len(Variants[v_id]['evidence_list'])
                    continue

                if Variants[v_id]['main_variant_class'] == MUTATION:
                    if Variants[v_id]['prot_ref_start_pos'] > 0:
                        if get_aachange_overlap_length( aachange, v_id, Variants ) > 0:
                            called_str = '{gene} {var}|{genomic_change}|{ref}'.format( gene=gene, var=aachange, genomic_change=gdnachange, ref=ref_build )
                            list_append( Matches[ sample_pair ]['partial'], {'v_id': v_id, 'reason': '4.pcoordOverlap', 'called': called_str} )
                            Variant_tracking[sample_pair][alteration_summary]['v_id_list'].append( v_id )
                            Variant_tracking[sample_pair][alteration_summary]['total_evidence_count'] += len(Variants[v_id]['evidence_list'])
                            continue

        # Check for potentially relevant trials
        if len(args.annotate_trials) and bHasSampleMatch:

            if re.search(r'[SDOTM]NP', vartype):
                if gene in Trials.keys():
                    if len( Trials[ gene ][ MUTATION ].keys() ):
                        GenesSeenInTrials[ MUTATION ].append( gene )

                        # check for nonspecific match
                        if 'any' in Trials[ gene ][ MUTATION ].keys():
                            check_alloc_named( Matches_trials[ sample_pair ][ MUTATION ], gene, 'list' )
                            Matches_trials[ sample_pair ][ MUTATION ][ gene ].extend( Trials[ gene ][ MUTATION ][ 'any' ] )
                            if myglobal.DEBUG:
                                logger.info('{gene} matched -any- trials for xNP' . format(gene=gene))

                        # check for aachanges "as is" or with 'any' alternative allele tag
                        p = aachange
                        m = re.search(r'([A-Z]\d+)(.*)', p )
                        any_probe = m.groups()[0] + ':any:'

                        if myglobal.DEBUG:
                            logger.info('we will check for vartype={x1}, gene={x2}, pos={x3}, any={x4}' . format(x1=vartype, x2=gene, x3=p, x4=any_probe))

                        if  p         in Trials[ gene ][ MUTATION ].keys():
                            check_alloc_named( Matches_trials[ sample_pair ][ MUTATION ], gene, 'list' )
                            Matches_trials[ sample_pair ][ MUTATION ][ gene ].extend( Trials[ gene ][ MUTATION ][ p ] )
                            if myglobal.DEBUG:
                                logger.info('...matched aachange as {}' . format(p))

                        if any_probe in Trials[ gene ][ MUTATION ].keys():
                            check_alloc_named( Matches_trials[ sample_pair ][ MUTATION ], gene, 'list' )
                            Matches_trials[ sample_pair ][ MUTATION ][ gene ].extend( Trials[ gene ][ MUTATION ][ any_probe ] )
                            if myglobal.DEBUG:
                                logger.info('...matched aachange as any probe {}' . format(any_probe))

                        # check for genomic coordinate
                        gc_probe = ':'.join([ chrom, pos_start, pos_end, ref, alt ])
                        if gc_probe in Trials[ gene ][ MUTATION ].keys():
                            check_alloc_named( Matches_trials[ sample_pair ][ MUTATION ], gene, 'list' )
                            Matches_trials[ sample_pair ][ MUTATION ][ gene ].extend( Trials[ gene ][ MUTATION ][ gc_probe ] )
                            if myglobal.DEBUG:
                                logger.info('...matched aachange as any probe {}' . format(gc_probe))


            elif vartype in ['INS', 'DEL']:
                for vt in [ INSERTION, DELETION ]:
                    if gene in Trials.keys():
                        if len( Trials[ gene ][ vt ].keys() ):
                            GenesSeenInTrials[ vt ].append( gene )

                            # check for nonspecific match
                            if 'any' in Trials[ gene ][ vt ].keys():
                                check_alloc_named( Matches_trials[ sample_pair ][ vt ], gene, 'list' )
                                Matches_trials[ sample_pair ][ vt ][ gene ].extend( Trials[ gene ][ vt ][ 'any' ] )
                                if myglobal.DEBUG:
                                    logger.info('{gene} matched -any- trials for {vartype}' . format(gene=gene, vartype=map_mut_reverse(vt)))

                            # check for aachanges "as is"
                            if vt == INSERTION:
                                if alt in Trials[ gene ][ vt ].keys():
                                    check_alloc_named( Matches_trials[ sample_pair ][ vt ], gene, 'list' )
                                    Matches_trials[ sample_pair ][ vt ][ gene ].extend( Trials[ gene ][ vt ][ alt ] )
                                    if myglobal.DEBUG:
                                        logger.info('{gene} matched trials for {vartype} with alteration {alteration}' . format(gene=gene, vartype=map_mut_reverse(vt), alteration=alt))

                            elif vt == DELETION:
                                if ref in Trials[ gene ][ vt ].keys():
                                    check_alloc_named( Matches_trials[ sample_pair ][ vt ], gene, 'list' )
                                    Matches_trials[ sample_pair ][ vt ][ gene ].extend( Trials[ gene ][ vt ][ ref ] )
                                    if myglobal.DEBUG:
                                        logger.info('{gene} matched trials for {vartype} with alteration {alteration}' . format(gene=gene, vartype=map_mut_reverse(vt), alteration=ref))

                            else:
                                pass

            else:
                logger.warning('Unexpected vartype {} in input variant file' . format(vartype) )

    tsv_file.close()

    # Record sample status
    if bHasSampleMatch:
        logger.info('Sample name was mentioned in the input file')
    else:
        logger.info('Sample name was NOT mentioned in the input file')

    # Other annotations of trials
    if len(args.annotate_trials) and bHasSampleMatch:
        Sample = sample_pair

        # Wildtypes
        evaluate_trials_wildtype( Trials, Genes_altered, Sample, Matches_trials, 'somatic' )

        # Evaluate gene lists related to trials
        GenesSeenInTrials['all_types'] = []    # one-off key to store combined gene list
        for var_type in [ MUTATION, INSERTION, DELETION ]:
            GenesSeenInTrials[ var_type ] = uniquify(GenesSeenInTrials[ var_type ])
            GenesSeenInTrials['all_types'].extend( GenesSeenInTrials[ var_type ] )
        GenesSeenInTrials['all_types'] = uniquify(GenesSeenInTrials['all_types'])
        mut_seen  =  ','.join(sorted( GenesSeenInTrials[MUTATION]  ))
        ins_seen  =  ','.join(sorted( GenesSeenInTrials[INSERTION] ))
        del_seen  =  ','.join(sorted( GenesSeenInTrials[DELETION]  ))

        wt_genes_trials = []
        for xgene in Trials.keys():
            if len( Trials[ xgene ][ WILDTYPE ].keys() ):
                wt_genes_trials.append( xgene )

        wt_seen         = ','.join(sorted( filter(lambda i: i not in GenesSeenInTrials['all_types'], wt_genes_trials) ))
        logger.info('Altered genes having MUTATIONS allowed by trials:'  + (mut_seen if len(mut_seen) else 'n/a'))
        logger.info('Altered genes having INSERTIONS allowed by trials:' + (ins_seen if len(ins_seen) else 'n/a'))
        logger.info('Altered genes having DELETIONS allowed by trials:'  + (del_seen if len(del_seen) else 'n/a'))
        logger.info('Non-altered genes allowed by trials:'               + (wt_seen  if len(wt_seen ) else 'n/a'))

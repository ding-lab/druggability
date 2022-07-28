
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


def process_maf( args, Matches, Evidence, Variants, Genes, Fasta, Genes_altered, Trials, Matches_trials ):
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
    for alt_type in [ MUTATION, INDEL ]:
        GenesSeenInTrials[ alt_type ] = []

    if len(args.annotate_trials):
        check_alloc_named( Matches_trials, sample_pair, 'dict')
        for vc in VARIANT_CLASSES:
            check_alloc_named( Matches_trials[ sample_pair ], vc, 'dict' )

    # Main
    for row in read_tsv:
        fields = [ s.strip() for s in row ]

        if bReadHeader:
            if row[0] != 'Hugo_Symbol':
                continue
            if ['Chromosome', 'Start_Position', 'End_Position', 'Variant_Type', 'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode', 'HGVSp_Short'] == [ fields[x] for x in [4,5,6,9,15,16,36] ]:
                maf_filetype = WASHU_MAF
                bReadHeader = False
                continue
            elif ['Chromosome', 'Start_Position', 'End_Position', 'Variant_Type', 'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode', 'Protein_Change'] == [ fields[x] for x in [3,4,5,8,12,13,26] ]:
                maf_filetype = UNION_MAF
                bReadHeader = False
                continue
            else:
                abort_run('Unrecognized maf format')


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
                    check_alloc_named( Matches, sample_pair, 'match_level' )
                    called_str = '{gene} {var}|{genomic_change}|{ref}'.format( gene=gene, var=aachange, genomic_change=gdnachange, ref=ref_build )
                    list_append( Matches[ sample_pair ]['full'], {'v_id': v_id, 'reason': '1.gcoordExact_altExact', 'called': called_str} )
                    Variant_tracking[sample_pair][alteration_summary]['v_id_list'].append( v_id )
                    Variant_tracking[sample_pair][alteration_summary]['total_evidence_count'] += len(Variants[v_id]['evidence_list'])
                    continue

                if has_genomic_match( gdnacoords, v_id, Variants, 'gdnacoords_liftover' ):
                    check_alloc_named( Matches, sample_pair, 'match_level' )
                    called_str = '{gene} {var}|{genomic_change}|{ref}'.format( gene=gene, var=aachange, genomic_change=gdnachange, ref=ref_build )
                    list_append( Matches[ sample_pair ]['partial'], {'v_id': v_id, 'reason': '2.gcoordExact_altDiff', 'called': called_str} )
                    Variant_tracking[sample_pair][alteration_summary]['v_id_list'].append( v_id )
                    Variant_tracking[sample_pair][alteration_summary]['total_evidence_count'] += len(Variants[v_id]['evidence_list'])
                    continue

                if has_genomic_overlap( tmp_set, v_id, Variants ):
                    check_alloc_named( Matches, sample_pair, 'match_level' )
                    called_str = '{gene} {var}|{genomic_change}|{ref}'.format( gene=gene, var=aachange, genomic_change=gdnachange, ref=ref_build )
                    list_append( Matches[ sample_pair ]['partial'], {'v_id': v_id, 'reason': '3.gcoordOverlap', 'called': called_str} )
                    Variant_tracking[sample_pair][alteration_summary]['v_id_list'].append( v_id )
                    Variant_tracking[sample_pair][alteration_summary]['total_evidence_count'] += len(Variants[v_id]['evidence_list'])
                    continue

                if Variants[v_id]['main_variant_class'] == MUTATION:
                    if Variants[v_id]['prot_ref_start_pos'] > 0:
                        if get_aachange_overlap_length( aachange, v_id, Variants ) > 0:
                            check_alloc_named( Matches, sample_pair, 'match_level' )
                            called_str = '{gene} {var}|{genomic_change}|{ref}'.format( gene=gene, var=aachange, genomic_change=gdnachange, ref=ref_build )
                            list_append( Matches[ sample_pair ]['partial'], {'v_id': v_id, 'reason': '4.pcoordOverlap', 'called': called_str} )
                            Variant_tracking[sample_pair][alteration_summary]['v_id_list'].append( v_id )
                            Variant_tracking[sample_pair][alteration_summary]['total_evidence_count'] += len(Variants[v_id]['evidence_list'])
                            continue

        # Check for potentially relevant trials
        if len(args.annotate_trials) and bHasSampleMatch:

            if( re.search(r'[SDOTM]NP', vartype) ):
                if( gene in Trials[ MUTATION ].keys() ):
                    GenesSeenInTrials[ MUTATION ].append( gene )

                    # check for 'any'
                    if( 'any' in Trials[ MUTATION ][ gene ].keys() ):
                        check_alloc_named( Matches_trials[ sample_pair ][ MUTATION ], gene, 'list' )
                        Matches_trials[ sample_pair ][ MUTATION ][ gene ].extend( Trials[ MUTATION ][ gene ][ 'any' ] )
                        if myglobal.DEBUG:
                            logger.info('{gene} matched -any- trials for xNP' . format(gene=gene))

                    # also check for aachanges "as is" or with 'any' alternative allele
                    for p in clean_split( aachange ):
                        m = re.search(r'([A-Z]\d+)(.*)', p )
                        any_probe = m.groups()[0] + ':any:'

                        if myglobal.DEBUG:
                            logger.info('we will check for vartype={x1}, gene={x2}, pos={x3}, any={x4}' . format(x1=vartype, x2=gene, x3=p, x4=any_probe))

                        if( p         in Trials[ MUTATION ][ gene ].keys() ):
                            check_alloc_named( Matches_trials[ sample_pair ][ MUTATION ], gene, 'list' )
                            Matches_trials[ sample_pair ][ MUTATION ][ gene ].extend( Trials[ MUTATION ][ gene ][ p ] )
                            if myglobal.DEBUG:
                                logger.info('...matched aachange as {}' . format(p))

                        if( any_probe in Trials[ MUTATION ][ gene ].keys() ):
                            check_alloc_named( Matches_trials[ sample_pair ][ MUTATION ], gene, 'list' )
                            Matches_trials[ sample_pair ][ MUTATION ][ gene ].extend( Trials[ MUTATION ][ gene ][ any_probe ] )
                            if myglobal.DEBUG:
                                logger.info('...matched aachange as any probe {}' . format(any_probe))

            elif( vartype in ['INS', 'DEL'] ):
                if( gene in Trials[ INDEL ].keys() ):
                    GenesSeenInTrials[ INDEL ].append( gene )

                    # check for 'any'
                    if( 'any' in Trials[ INDEL ][ gene ].keys() ):
                        check_alloc_named( Matches_trials[ sample_pair ][ INDEL ], gene, 'list' )
                        Matches_trials[ sample_pair ][ INDEL ][ gene ].extend( Trials[ INDEL ][ gene ][ 'any' ] )
                        if myglobal.DEBUG:
                            logger.info('{gene} matched -any- trials for indel' . format(gene=gene))

                    # also check for aachanges "as is" or with 'any' alternative allele
                    # currently there are none such
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
        evaluate_trials_wildtype( Trials, Genes_altered, Sample, Matches_trials )

        # Evaluate gene lists related to trials
        GenesSeenInTrials['all_types'] = []    # one-off key
        for var_type in [ MUTATION, INDEL ]:
            GenesSeenInTrials[ var_type ] = uniquify(GenesSeenInTrials[ var_type ])
            GenesSeenInTrials['all_types'].extend( GenesSeenInTrials[ var_type ] )
        GenesSeenInTrials['all_types'] = uniquify(GenesSeenInTrials['all_types'])
        mut_seen        = ','.join(sorted( GenesSeenInTrials[MUTATION] ))
        indel_seen      = ','.join(sorted( GenesSeenInTrials[INDEL] ))
        wt_genes_trials = [ x for x in Trials[ WILDTYPE ].keys()]
        wt_seen         = ','.join(sorted( filter(lambda i: i not in GenesSeenInTrials['all_types'], wt_genes_trials) ))
        logger.info('Altered genes having MUTATIONS allowed by trials:' + (mut_seen   if len(mut_seen)   else 'n/a'))
        logger.info('Altered genes having INDELS allowed by trials:'    + (indel_seen if len(indel_seen) else 'n/a'))
        logger.info('Non-altered genes allowed by trials:'              + (wt_seen    if len(wt_seen)    else 'n/a'))

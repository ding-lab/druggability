
# R. Jay Mashl <rmashl@wustl.edu>

import os, sys, csv, re
import config
from utils import *
from enums import *
from harmonize import *
import logging

DEBUG=config.DEBUG
DEBUG_2=config.DEBUG_2

# accommodate large strings
csv.field_size_limit( 131072 * 82 )


def process_maf( args, Evidence, Variants, Genes, Fasta):

    inputFile = args.variant_file
    Variant_tracking = dict()   # record which samples have which variants
    Matches          = dict()   # matches by sample, separated in 'full' and 'partial' match lists
    maf_filetype     = UNDECLARED

    bReadHeader       = True

    tsv_file = open( inputFile )
    read_tsv = csv.reader(tsv_file, delimiter='\t')
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
                logging.error('Unrecognized maf format')
                sys.exit(1)

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
            vartype      = fields[ 8]   # SNP, INS, DEL, ...

            ref = fields[ 9]

            if ref == fields[10]:
                alt = fields[11] if fields[11] != '-' else ''
            elif ref == fields[11]:
                alt = fields[10] if fields[10] != '-' else ''
            else:
                logging.error('cannot resolve alteration at {gene} {pos}'.format( gene=gene, pos=str(pos_start) ))
                sys.exit(1)

            sample_t     = fields[12]   # here, this is tumor sample
            sample_n     = fields[13]   # here, this is (matched) normal sample
            cdnachange   = fields[24]
            aachange     = fields[26]



        # Check whether this gene is mentioned in any database
        if gene not in Genes.keys():
            continue

        # summarize alteration
        alteration_summary = '\t'.join([ gene, chrom, pos_start, cdnachange, aachange, ref_build, maf_varclass, vartype ])

        # initially look at only vars with AA change in HGVS short format; if blank, it is often a splice site
        if not re.search( r'^p\.', aachange):
            if DEBUG_2:
                logging.info( 'maf record IGNORED: ' + alteration_summary )
            continue

        # harmonize
        aachange = harmonize_maf( aachange )

        if maf_filetype == UNION_MAF:       # additional treatment
            if re.search(r'>', aachange) or re.search('ins', aachange) or re.search(r'del$', aachange):
                aachange = harmonize_maf_2( aachange, gene, Fasta )

        # calculate gdna change
        tmp_set = { 'chrom': chrom, 'pos0': pos_start, 'pos1': pos_end, 'ref': ref, 'alt': alt, 'gene': gene }
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

        # Set up storage for tracking matches using composite key
        sample_pair = '{tumor}||{normal}'.format( tumor=sample_t, normal=sample_n )
        if sample_pair not in Variant_tracking.keys():
            Variant_tracking[ sample_pair ] = dict()
        Variant_tracking[ sample_pair ][alteration_summary] = dict( total_evidence_count=0, v_id_list=[] )

        # Identify and track matches
        num_full_matches    = 0
        num_partial_matches = 0
        num_unmatched       = 0


        for v_id in Genes[gene]:

            if has_genomic_match( gdnachange, v_id, Variants, 'gdnachange_liftover' ):
                check_alloc_match( Matches, sample_pair )
                called_str = '{gene} {var}|{genomic_change}|{ref}'.format( gene=gene, var=aachange, genomic_change=gdnachange, ref=ref_build )
                list_append( Matches[ sample_pair ]['full'], {'v_id': v_id, 'reason': '1.gcoordExact_altExact', 'called': called_str} )
                Variant_tracking[sample_pair][alteration_summary]['v_id_list'].append( v_id )
                Variant_tracking[sample_pair][alteration_summary]['total_evidence_count'] += len(Variants[v_id]['evidence_list'])
                continue

            if has_genomic_match( gdnacoords, v_id, Variants, 'gdnacoords_liftover' ):
                check_alloc_match( Matches, sample_pair )
                called_str = '{gene} {var}|{genomic_change}|{ref}'.format( gene=gene, var=aachange, genomic_change=gdnachange, ref=ref_build )
                list_append( Matches[ sample_pair ]['partial'], {'v_id': v_id, 'reason': '2.gcoordExact_altDiff', 'called': called_str} )
                Variant_tracking[sample_pair][alteration_summary]['v_id_list'].append( v_id )
                Variant_tracking[sample_pair][alteration_summary]['total_evidence_count'] += len(Variants[v_id]['evidence_list'])
                continue

            if has_genomic_overlap( tmp_set, v_id, Variants ):
                check_alloc_match( Matches, sample_pair )
                called_str = '{gene} {var}|{genomic_change}|{ref}'.format( gene=gene, var=aachange, genomic_change=gdnachange, ref=ref_build )
                list_append( Matches[ sample_pair ]['partial'], {'v_id': v_id, 'reason': '3.gcoordOverlap', 'called': called_str} )
                Variant_tracking[sample_pair][alteration_summary]['v_id_list'].append( v_id )
                Variant_tracking[sample_pair][alteration_summary]['total_evidence_count'] += len(Variants[v_id]['evidence_list'])
                continue

            if Variants[v_id]['main_variant_class'] == MUTATION:
                if Variants[v_id]['prot_ref_start_pos'] > 0:
                    if get_aachange_overlap_length( aachange, v_id, Variants ) > 0:
                        check_alloc_match( Matches, sample_pair )
                        called_str = '{gene} {var}|{genomic_change}|{ref}'.format( gene=gene, var=aachange, genomic_change=gdnachange, ref=ref_build )
                        list_append( Matches[ sample_pair ]['partial'], {'v_id': v_id, 'reason': '4.pcoordOverlap', 'called': called_str} )
                        Variant_tracking[sample_pair][alteration_summary]['v_id_list'].append( v_id )
                        Variant_tracking[sample_pair][alteration_summary]['total_evidence_count'] += len(Variants[v_id]['evidence_list'])
                        continue



    tsv_file.close()


    # Print summary
    #print_summary_by_sample( Variant_tracking, Variants, Evidence )

    print_summary_for_all( Matches, Variants, Evidence, args )

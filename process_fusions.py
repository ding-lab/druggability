
# R. Jay Mashl <rmashl@wustl.edu>

import os, sys, csv, re
import druggability_databases.config as config
import myglobal
from utils import *
from enums import *
import logging

def process_fusions( args, Evidence, Variants, Genes):
    inputFile         = args.variant_file
    Variant_tracking  = dict()   # record variants by sample
    fusion_filetype   = UNDECLARED
    Matches           = dict()   # matches by sample, separated in 'full' and 'partial' match lists

    hdr               = []    # column headers
    expected_hdr      = ['FusionName', 'LeftBreakpoint', 'RightBreakpoint', 'Sample', 'JunctionReadCount', 'SpanningFragCount', 'FFPM', 'PROT_FUSION_TYPE']
    bReadHeader       = True
    bHasSampleMatch   = False

    tsv_file  = open( inputFile )
    read_tsv  = csv.reader(tsv_file, delimiter='\t')
    for row in read_tsv:
        fields = [ s.strip() for s in row ]

        if bReadHeader:
            hdr = fields
            if row[0] == 'FusionName':
                fusion_filetype = SINGLE
                bReadHeader = False
                if hdr[0:8] != expected_hdr:
                    abort_run('Unexpected fusion file column format')
                continue
            elif row[2] == 'FusionName':
                fusion_filetype = COMBINED
                bReadHeader = False
                if hdr[2:10] != expected_hdr:
                    abort_run('Unexpected fusion file column format')
                continue
            else:
                abort_run('Unrecognized fusion file format')

        # shift entries for aggregated fusion reports to match single-sample reports
        if fusion_filetype == COMBINED:
            fields = fields[2:]

        # Process entry

        FusionName, LeftBreakpoint, RightBreakpoint, Sample, JunctionReadCount, SpanningFragCount, FFPM, PROT_FUSION_TYPE = fields[0:8]

        # Require sample to match
        if Sample != args.tumor_name:
            continue
        bHasSampleMatch = True

        # summarize alteration
        alteration_summary = '\t'.join([ FusionName, LeftBreakpoint, RightBreakpoint ])

        # Set up storage for tracking matches
        if Sample not in Variant_tracking.keys():
            Variant_tracking[ Sample ] = dict()
        Variant_tracking[ Sample ][alteration_summary] = dict( total_evidence_count=0, v_id_list=[] )

        # Gather gene names
        genes = []
        for g in FusionName.split('--'):
            genes.append(  g.split('.')[0] )   # remove transcript identifier

        # Identify matches
        for g in genes:
            if g in Genes.keys():
                for v_id in Genes[g]:
                    if  Variants[v_id]['main_variant_class'] == FUSION:

                        # Determine full vs partial gene match (wildcard gene is handled implicitly)
                        myoverlap = intersection( genes, Variants[v_id]['fusion_gene_set'] )
                        num_hits = len( intersection( genes, Variants[v_id]['fusion_gene_set'] ) )
                        if num_hits == 2:

                            if Variants[v_id]['pp_conditions'] == 1:   # no additional criteria
                                check_alloc_match( Matches, Sample )
                                list_append( Matches[ Sample ]['full'], {'v_id': v_id, 'reason': '-', 'called': FusionName} )   # unchecked criteria
                            else:
                                check_alloc_match( Matches, Sample )
                                list_append( Matches[ Sample ]['partial'], {'v_id': v_id, 'reason': Variants[v_id]['pp_condition2_value'], 'called': FusionName} )

                        elif num_hits == 1:
                            check_alloc_match( Matches, Sample )
                            list_append( Matches[ Sample ]['partial'], {'v_id': v_id, 'reason': 'matched ' + str(myoverlap[0]), 'called': FusionName} )
                        else:
                            pass


                        Variant_tracking[Sample][alteration_summary]['v_id_list'].append( v_id )
                        Variant_tracking[Sample][alteration_summary]['total_evidence_count'] += len(Variants[v_id]['evidence_list'])

    tsv_file.close()

    if bHasSampleMatch:
        logging.info('Sample name was mentioned in the input file')
    else:
        logging.info('Sample name was NOT mentioned in the input file')

    # Print summary by sample
    # print_summary_by_sample( Variant_tracking, Variants, Evidence )

    print_summary_for_all( Matches, Variants, Evidence, args )

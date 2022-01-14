
# R. Jay Mashl <rmashl@wustl.edu>

import os, sys, csv, re
import config
from utils import *

DEBUG=config.DEBUG

def process_fusions( args, Evidence, Variants, Genes):

    inputFile = args.variant_file
    tsv_file = open( inputFile )
    read_tsv = csv.reader(tsv_file, delimiter='\t')

    gene_evidence = dict()       # format:  gene => evidence[]
    Variant_sample_tracking = dict()   # record which samples have which variants

    # fusion filetypes
    SINGLE   = 1
    COMBINED = 2
    fusion_filetype = 0

    col = []    # field names
    bReadHeader = True
    for row in read_tsv:
        fields = [ s.strip() for s in row ]

        if bReadHeader:
            if row[0] == "FusionName":
                fusion_filetype = SINGLE
                bReadHeader = False
                continue
            elif row[2] == "FusionName":
                fusion_filetype = COMBINED
                bReadHeader = False
                continue
            else:
                print("ERROR: Unrecognized fusion file format")
                sys.exit(1)

        # shift entries for aggregated fusion reports and leave single-sample reports alone
        if fusion_filetype == COMBINED:
            fields = fields[2:]

        FusionName, LeftBreakpoint, RightBreakpoint, Sample, JunctionReadCount, SpanningFragCount, FFPM, PROT_FUSION_TYPE = fields[0:8]

        # summarize alteration
        alteration_summary = '\t'.join([ FusionName, LeftBreakpoint, RightBreakpoint ])

        # Set up storage for tracking matches
        if Sample not in Variant_sample_tracking.keys():
            Variant_sample_tracking[ Sample ] = dict()
        Variant_sample_tracking[ Sample ][alteration_summary] = dict( total_evidence_count=0, v_id_list=[] )

        # Identify matches
        genes = []
        for g in FusionName.split('--'):
            genes.append(  g.split('.')[0] )   # remove transcript identifier

        for g in genes:
            if g in Genes.keys():
                for v_id in Genes[g]:
                    if  Variants[v_id]['our_variant_category'] == 'fusion':
                        Variant_sample_tracking[Sample][alteration_summary]['v_id_list'].append( v_id )
                        Variant_sample_tracking[Sample][alteration_summary]['total_evidence_count'] += len(Variants[v_id]['evidence_list'])

    tsv_file.close()


    # Print summary by sample
    for sample in Variant_sample_tracking.keys():
        for alteration in Variant_sample_tracking[sample].keys():
            this_alt = Variant_sample_tracking[sample][alteration]
            if this_alt['total_evidence_count']:
                print_sample_header( sample, alteration )
                print_output_header()
                for v_id in this_alt['v_id_list']:
                    for ev_id in Variants[v_id]['evidence_list']:
                        t = Evidence[ev_id]
                        print( *[ (v_id.split(':'))[0],   Variants[v_id]['variant'], t['disease'], t['oncogenicity'], t['mutation_effect'],   t['drugs_list_string'], t['evidence_type'], t['evidence_direction'], t['evidence_level'], t['clinical_significance'], format_citations(t['citations'])], sep = '\t')

                print('')
                print('')


# R. Jay Mashl <rmashl@wustl.edu>

import csv, re
import druggability_databases.config as config
import myglobal
from utils import *
from enums import *
import logging
from annotate_trials import *

logger = logging.getLogger(__name__)
logger.setLevel(0)

def process_fusions( args, Matches, Evidence, Variants, Genes, Genes_altered, Trials, Matches_trials, SampleMentioned, call_context, GenesSeenInTrials ):
    inputFile         = args.variant_fusion_file
    #Variant_tracking  = dict()   # record variants by sample
    fusion_filetype   = UNDECLARED

    hdr               = []    # column headers
    expected_hdr      = ['FusionName', 'LeftBreakpoint', 'RightBreakpoint', 'Sample', 'JunctionReadCount', 'SpanningFragCount', 'FFPM', 'PROT_FUSION_TYPE']
    bReadHeader       = True
    #bHasSampleMatch   = False
    SampleMentioned['fusion'] = False

    tsv_file  = open( inputFile )
    read_tsv  = csv.reader(tsv_file, delimiter='\t')

    # set key for tracking matches
    #Sample = args.tumor_name
    Sample = SAMPLENAME


    # Main
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
        FusionName, LeftBreakpoint, RightBreakpoint, this_Sample, JunctionReadCount, SpanningFragCount, FFPM, PROT_FUSION_TYPE = fields[0:8]

        # Require sample to match
        if this_Sample != args.fusion_sample_name:
            continue
        #bHasSampleMatch = True
        SampleMentioned['fusion'] = True

        # summarize alteration
        alteration_summary = '\t'.join([ FusionName, LeftBreakpoint, RightBreakpoint ])
        called_str = 'FUS:' + FusionName

        # Set up storage for tracking matches with alteration db
        check_alloc_named( Matches, Sample, 'match_level' )

        # Gather gene names
        genes = []
        for g in FusionName.split('--'):
            genes.append(  g.split('.')[0] )   # remove transcript identifier

        # Identify matches with alteration db
        for g in genes:
            Genes_altered[ g ] = 1
            if g in Genes.keys():
                for v_id in Genes[g]:
                    if  Variants[v_id]['main_variant_class'] == FUSION:

                        # Determine full vs partial gene match (wildcard gene is handled implicitly)
                        # The field 'reason' corresponds to criteria met and currently refers to gene-level matches
                        myoverlap = intersection( genes, Variants[v_id]['fusion_gene_set'] )
                        num_hits  = len( myoverlap )
                        if num_hits == 2:
                            if Variants[v_id]['pp_conditions'] == 1:   # no additional criteria
                                list_append( Matches[ Sample ]['full'], {'v_id': v_id, 'reason': '1.druggable gene pair::no additional criteria', 'called': called_str} )
                            else:
                                list_append( Matches[ Sample ]['partial'], {'v_id': v_id, 'reason': '2.druggable gene pair::additional criteria not applied', 'called': called_str} )

                        elif num_hits == 1:
                            if '*' in Variants[v_id]['fusion_gene_set']:  # wildcard present
                                if Variants[v_id]['pp_conditions'] == 1:   # no additional criteria
                                    list_append( Matches[ Sample ]['partial'], {'v_id': v_id, 'reason': '3.druggable gene with nonspecific partner::no additional criteria', 'called': called_str} )
                                else:
                                    list_append( Matches[ Sample ]['partial'], {'v_id': v_id, 'reason': '4.druggable gene with nonspecific partner::additional criteria not applied', 'called': called_str} )
                            else:  # just one partner matched
                                if Variants[v_id]['pp_conditions'] == 1:   # no additional criteria
                                    list_append( Matches[ Sample ]['partial'], {'v_id': v_id, 'reason': '5.possibly druggable gene with nonmatching partner::no additional criteria', 'called': called_str} )
                                else:
                                    list_append( Matches[ Sample ]['partial'], {'v_id': v_id, 'reason': '6.possibly druggable gene with nonmatching partner::additional criteria not applied', 'called': called_str} )

                        else:
                            pass

            # Check for potentially relevant trials
            # Fusion search is simpler because only a gene-level match is needed; do the search below instead
    tsv_file.close()

    if bReadHeader:
        abort_run('Variant file {} is missing the header line' . format( inputFile ))

    # Evaluate gene lists related to trials
    if args.annotate_trials:
        for gene in Genes_altered.keys():
            if gene in Trials.keys():
                pos = 'any'
                if call_context in Trials[ gene ][ FUSION ].keys():
                    if len( Trials[ gene ][ FUSION ][ call_context ].keys() ):
                        GenesSeenInTrials[ FUSION ].append( gene )
                        if pos in Trials[ gene ][ FUSION ][ call_context ].keys():
                            for t in Trials[ gene ][ FUSION ][ call_context ][ pos ]:
                                check_alloc_named( Matches_trials, Sample, 'dict' )
                                check_alloc_named( Matches_trials[ Sample ], FUSION, 'dict' )
                                check_alloc_named( Matches_trials[ Sample ][ FUSION ], gene, 'list' )
                                Matches_trials[ Sample ][ FUSION ][ gene ].append( t )

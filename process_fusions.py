
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

def process_fusions( args, Matches, Evidence, Variants, Genes, Genes_altered, Trials, Matches_trials, call_mode ):
    inputFile         = args.variant_file
    Variant_tracking  = dict()   # record variants by sample
    fusion_filetype   = UNDECLARED

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

        # Set up storage for tracking matches with alteration db
        check_alloc_named( Matches, Sample, 'match_level' )

        if Sample not in Variant_tracking.keys():
            Variant_tracking[ Sample ] = dict()
        Variant_tracking[ Sample ][alteration_summary] = dict( total_evidence_count=0, v_id_list=[] )

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
                                list_append( Matches[ Sample ]['full'], {'v_id': v_id, 'reason': '1.druggable gene pair::no additional criteria', 'called': FusionName} )
                            else:
                                list_append( Matches[ Sample ]['partial'], {'v_id': v_id, 'reason': '2.druggable gene pair::additional criteria not applied', 'called': FusionName} )

                        elif num_hits == 1:
                            if '*' in Variants[v_id]['fusion_gene_set']:  # wildcard present
                                if Variants[v_id]['pp_conditions'] == 1:   # no additional criteria
                                    list_append( Matches[ Sample ]['partial'], {'v_id': v_id, 'reason': '3.druggable gene with nonspecific partner::no additional criteria', 'called': FusionName} )
                                else:
                                    list_append( Matches[ Sample ]['partial'], {'v_id': v_id, 'reason': '4.druggable gene with nonspecific partner::additional criteria not applied', 'called': FusionName} )
                            else:  # just one partner matched
                                if Variants[v_id]['pp_conditions'] == 1:   # no additional criteria
                                    list_append( Matches[ Sample ]['partial'], {'v_id': v_id, 'reason': '5.possibly druggable gene with nonmatching partner::no additional criteria', 'called': FusionName} )
                                else:
                                    list_append( Matches[ Sample ]['partial'], {'v_id': v_id, 'reason': '6.possibly druggable gene with nonmatching partner::additional criteria not applied', 'called': FusionName} )

                        else:
                            pass

                        Variant_tracking[Sample][alteration_summary]['v_id_list'].append( v_id )
                        Variant_tracking[Sample][alteration_summary]['total_evidence_count'] += len(Variants[v_id]['evidence_list'])

    tsv_file.close()

    # Record sample status
    if bHasSampleMatch:
        logger.info('Sample name was mentioned in the input file')
    else:
        logger.info('Sample name was NOT mentioned in the input file')

    # Evaluate gene lists related to trials
    if len(args.annotate_trials) and bHasSampleMatch:
        Sample = args.tumor_name

        # Set up storage
        check_alloc_named( Matches_trials, Sample, 'dict')
        for vc in VARIANT_CLASSES:
            check_alloc_named( Matches_trials[ Sample ], vc, 'dict' )

        evaluate_trials_fusion( Trials, Genes_altered, Sample, Matches_trials, call_mode )
        evaluate_trials_wildtype( Trials, Genes_altered, Sample, Matches_trials, call_mode )
        fusion_seen = ','.join(sorted(Matches_trials[ Sample ][ FUSION   ].keys()))
        wt_seen     = ','.join(sorted(Matches_trials[ Sample ][ WILDTYPE ].keys()))
        logger.info('Altered genes having FUSIONS allowed by trials:' + (fusion_seen if len(fusion_seen) else 'n/a'))
        logger.info('Non-altered genes allowed by trials:'            + (wt_seen     if len(wt_seen)     else 'n/a'))

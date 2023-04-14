
# R. Jay Mashl <rmashl@wustl.edu>

import os, sys, csv, re
import druggability_databases.config as config
import myglobal
from utils import *
from enums import *
import logging

logger = logging.getLogger(__name__)
logger.setLevel(0)

# Find trials mentioning gene wild types that we found as the reference allele
def evaluate_trials_wildtype( Trials, Genes_altered, Sample, Matches_trials, call_context, wt_seen ):
    for gene in Trials.keys():
        if gene not in Genes_altered.keys():  # wildtype genes in trials
            pos = 'none'
            if call_context in Trials[ gene ][ WILDTYPE ].keys():
                if pos in Trials[ gene ][ WILDTYPE ][ call_context ].keys():
                    wt_seen.append( gene )
                    for t in Trials[ gene ][ WILDTYPE ][ call_context ][ pos ]:
                        check_alloc_named( Matches_trials, Sample, 'dict' )
                        check_alloc_named( Matches_trials[ Sample ], WILDTYPE, 'dict' )
                        check_alloc_named( Matches_trials[ Sample ][ WILDTYPE ], gene, 'list' )
                        Matches_trials[ Sample ][ WILDTYPE ][ gene ].append( t )

# Find qualifications among trials
def evaluate_trials_disqualifications( Matches_trials ):
    for sample in Matches_trials.keys():
        Matches_trials[ sample ][ DISQUALIFYING] = []   # one-off key for listing disqualified trials
        for vt in VARIANT_CLASSES:
            if len(Matches_trials[ sample ][ vt ].keys()):
                for gene in Matches_trials[ sample ][ vt ].keys():
                    for hit in Matches_trials[ sample ][ vt ][ gene ]:
                        if hit['eligibility_type'] == DISQUALIFYING:
                            Matches_trials[ sample ][ DISQUALIFYING ].append( {'trial_id': hit['trial_id'], 'reason': hit['position_target']} )

# Create list of trials mentioning wildtype
def create_list_of_trials_with_wildtype( Trials, wt_genes_trials, call_context ):
    for xgene in Trials.keys():
        if call_context in Trials[ xgene ][ WILDTYPE ].keys():
            if len( Trials[ xgene ][ WILDTYPE ][ call_context ].keys() ):
                wt_genes_trials.append( xgene )

# Process trials
def evaluate_gene_trials_hits( args, Trials, Genes_altered, SAMPLENAME, Matches_trials, call_context, GenesSeenInTrials ):

    # Check for wildtypes
    wt_seen = []
    evaluate_trials_wildtype( Trials, Genes_altered, SAMPLENAME, Matches_trials, call_context, wt_seen)
    wt_seen = sorted(wt_seen)

    # Create single-type and merged list of altered genes in input(s) that are the subject of trials
    for var_type in [ MUTATION, INSERTION, DELETION, FUSION ]:
        GenesSeenInTrials[ var_type ] = uniquify(GenesSeenInTrials[ var_type ])
        GenesSeenInTrials['all_types'].extend( GenesSeenInTrials[ var_type ] )
    GenesSeenInTrials['all_types'] = uniquify(GenesSeenInTrials['all_types'])
    all_seen    = ','.join(sorted( GenesSeenInTrials['all_types'] ))
    mut_seen    = ','.join(sorted( GenesSeenInTrials[MUTATION] ))
    ins_seen    = ','.join(sorted( GenesSeenInTrials[INSERTION] ))
    del_seen    = ','.join(sorted( GenesSeenInTrials[DELETION] ))
    fusion_seen = ','.join(sorted( GenesSeenInTrials[FUSION] ))

    logger.info('Altered genes evaluated for trial matching in {} context: {}' . format(call_context, all_seen if len(all_seen) else 'n/a'))
    if args.variant_maf_file or args.variant_basicmaf_file:
        logger.info('Altered genes with MUTATIONS evaluated for trial matching in {} context: {}' . format(call_context, mut_seen if len(mut_seen) else 'n/a'))
        logger.info('Altered genes with INSERTIONS evaluated for trial matching in {} context: {}' . format(call_context, ins_seen if len(ins_seen) else 'n/a'))
        logger.info('Altered genes with DELETIONS evaluated for trial matching in {} context: {}' .format(call_context, del_seen if len(del_seen) else 'n/a'))

    if args.variant_fusion_file:
        logger.info('Altered genes with FUSIONS evaluated for trials matching in {} context: {}' . format(call_context, fusion_seen if len(fusion_seen) else 'n/a'))

    logger.info('Unaltered genes evaluated for trial matching in {} context: {}' . format(call_context, ','.join(wt_seen) if len(wt_seen) else 'n/a'))


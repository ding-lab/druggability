
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
                if pos in Trials[ gene ][ WILDTYPE ][ call_context ].keys():   # requires all WT-related info be mentioned in trials db
                    wt_seen.append( gene )
                    for t in Trials[ gene ][ WILDTYPE ][ call_context ][ pos ]:
                        check_alloc_named( Matches_trials, Sample, 'dict' )
                        check_alloc_named( Matches_trials[ Sample ], WILDTYPE, 'dict' )
                        check_alloc_named( Matches_trials[ Sample ][ WILDTYPE ], gene, 'list' )
                        Matches_trials[ Sample ][ WILDTYPE ][ gene ].append( t )

# Find trials mentioning a given gene name
def evaluate_trials_fusion( Trials, Genes_altered, Sample, Matches_trials, call_context ):
    for gene in Genes_altered.keys():
        if gene in Trials.keys():
            pos = 'any'
            if call_context in Trials[ gene ][ FUSION ].keys():
                if pos in Trials[ gene ][ FUSION ][ call_context ].keys():
                    for t in Trials[ gene ][ FUSION ][ call_context ][ pos ]:
                        check_alloc_named( Matches_trials, Sample, 'dict' )
                        check_alloc_named( Matches_trials[ Sample ], FUSION, 'dict' )
                        check_alloc_named( Matches_trials[ Sample ][ FUSION ], gene, 'list' )
                        Matches_trials[ Sample ][ FUSION ][ gene ].append( t )

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

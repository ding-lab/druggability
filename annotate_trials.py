
# R. Jay Mashl <rmashl@wustl.edu>

import os, sys, csv, re
import druggability_databases.config as config
import myglobal
from utils import *
from enums import *
import logging

logger = logging.getLogger(__name__)
logger.setLevel(0)

# Find trials mentioning gene wild types
def evaluate_trials_wildtype( Trials, Genes_altered, Sample, Matches_trials, call_mode ):
    for gene in Trials.keys():
        if gene not in Genes_altered.keys():  # wildtype genes in trials
            pos = 'none'
            if pos in Trials[ gene ][ WILDTYPE ].keys():
                for t in Trials[ gene ][ WILDTYPE ][ pos ]:
                    if t['call_mode'] == call_mode:
                        check_alloc_named( Matches_trials, Sample, 'dict' )
                        check_alloc_named( Matches_trials[ Sample ], WILDTYPE, 'dict' )
                        check_alloc_named( Matches_trials[ Sample ][ WILDTYPE ], gene, 'list' )
                        Matches_trials[ Sample ][ WILDTYPE ][ gene ].append( t )


# Find trials mentioning a given gene name
def evaluate_trials_fusion( Trials, Genes_altered, Sample, Matches_trials, call_mode ):
    for gene in Genes_altered.keys():
        if gene in Trials.keys():
            pos = 'any'
            if pos in Trials[ gene ][ FUSION ].keys():
                for t in Trials[ gene ][ FUSION ][ pos ]:
                    if t['call_mode'] == call_mode:
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

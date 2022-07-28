
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
def evaluate_trials_wildtype( Trials, Genes_altered, Sample, Matches_trials ):
    for wt_gene in Trials[ WILDTYPE ].keys():
        if wt_gene not in Genes_altered.keys():
            pos = 'none'
            for t in Trials[ WILDTYPE ][ wt_gene ][ pos ]:
                check_alloc_named( Matches_trials, Sample, 'dict' )
                check_alloc_named( Matches_trials[ Sample ], WILDTYPE, 'dict' )
                check_alloc_named( Matches_trials[ Sample ][ WILDTYPE ], wt_gene, 'list' )
                Matches_trials[ Sample ][ WILDTYPE ][ wt_gene ].append( t )


# Find trials mentioning a given gene name
def evaluate_trials_fusion( Trials, Genes_altered, Sample, Matches_trials ):
    for gene in Genes_altered.keys():
        if gene in Trials[ FUSION ].keys():
            pos = 'any'
            for t in Trials[ FUSION ][ gene ][ pos ]:
                check_alloc_named( Matches_trials, Sample, 'dict' )
                check_alloc_named( Matches_trials[ Sample ], FUSION, 'dict' )
                check_alloc_named( Matches_trials[ Sample ][ FUSION ], gene, 'list' )
                Matches_trials[ Sample ][ FUSION ][ gene ].append( t )

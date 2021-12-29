#!/usr/bin/env python3

# R. Jay Mashl <rmashl@wustl.edu>

import getopt
import os, sys, csv, re
import argparse

import config
from load_databases import *
from process_maf import *
from process_fusions import *

Evidence = dict()   # usage: key= drug, value = dict()
Variants       = dict()
Genes          = dict()   # index the variants
VariantAliases = dict()   # for possibly crosslinking databases

DEBUG=config.DEBUG


# Process command line
parser = argparse.ArgumentParser()
#parser.add_argument('-o', dest='output_file', type=str, required=False, help='output filename', default='stdout')
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('-t', dest='variation_type', type=str, required=True, help='variation type:  maf | fusion')
requiredNamed.add_argument('-f', dest='variant_file', type=str, required=True, help='variant filename')
args = parser.parse_args()
if args.variation_type not in ['maf', 'fusion']:
    print("ERROR: please specifiy variation type: maf | fusion")
    sys.exit(1)


# load drugbank
#drugbank_db = dict()
#with open( drugbank_files['file'], 'r') as json_file:
#    parsed = json.load( json_file )
#    for gene,value in parsed.items():
#        drugbank_db[ gene ] = value


# load variant summary from CIViC
load_civic( Variants, Genes, VariantAliases)
if DEBUG:
    print('# num variants loaded so far: ' + str(len(Variants.keys())))
    print('# num genes loaded so far: ' + str(len(Genes.keys())))

# load variant summary from OncoKB
load_oncokb( Variants, Genes, VariantAliases)
if DEBUG:
    print('# num variants loaded so far: ' + str(len(Variants.keys())))
    print('# num genes loaded so far: ' + str( len(Genes.keys())))

# Load evidence from CIViC
load_civic_evidence( Evidence, Variants )

# Load evidence1 from oncokb
load_oncokb_evidence( Evidence, Variants )
load_oncokb_therapeutics( Evidence, Variants, Genes )


if args.variation_type == 'maf':
    process_maf( args, Evidence, Variants, Genes)

if args.variation_type == 'fusion':
    process_fusions( args, Evidence, Variants, Genes)

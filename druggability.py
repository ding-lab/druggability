#!/usr/bin/env python3

# R. Jay Mashl <rmashl@wustl.edu>

import getopt
import os, sys, csv, re
import argparse
import logging

import config
from load_databases import *
from process_maf import *
from process_fusions import *
from utils import abort_run

def main( args ):
    Evidence       = dict()   # usage: key= drug, value = dict()
    Variants       = dict()
    Genes          = dict()   # index the variants
    VariantAliases = dict()   # for possibly crosslinking databases
    Fasta          = dict()   # fasta sequences for reformatting variants

    DEBUG=config.DEBUG

    # load variant summary from CIViC
    load_civic( Variants, Genes, VariantAliases)

    # Load evidence from CIViC
    load_civic_evidence( Evidence, Variants )

    if args.variation_type == 'maf':
        process_maf( args, Evidence, Variants, Genes, Fasta)

    if args.variation_type == 'fusion':
        process_fusions( args, Evidence, Variants, Genes)



if __name__ == '__main__':

    # Process command line
    parser = argparse.ArgumentParser()
    #parser.add_argument('-o', dest='output_file', type=str, required=False, help='output filename', default='stdout')
    parser.add_argument('-l', dest='log_file', type=str, required=False, help='logfile name', default='druggability.log')
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-t', dest='variation_type', type=str, required=True, help='variation type:  maf | fusion')
    requiredNamed.add_argument('-f', dest='variant_file', type=str, required=True, help='variant filename')
    args = parser.parse_args()

    # Set up log file
    mylogfile = args.log_file
    if os.path.exists( mylogfile ):
        os.remove( mylogfile )
    logging.basicConfig( filename = mylogfile, level=logging.DEBUG)

    # Validate input and run
    log_timestamp('start time')
    if args.variation_type not in ['maf', 'fusion']:
        abort_run('invalid input variation type')
    main( args )
    log_timestamp('end time')

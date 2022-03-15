#!/usr/bin/env python3

# R. Jay Mashl <rmashl@wustl.edu>

import os, sys
from pathlib import Path

# Update for module location
sys.path.append(Path(__file__).parent / "druggability_databases")

import getopt
import csv, re
import argparse
import logging

import myglobal
from load_databases import *
from process_maf import *
from process_fusions import *
from utils import abort_run
import druggability_databases.config as config

def main( args ):
    Evidence       = dict()   # usage: key= drug, value = dict()
    Variants       = dict()
    Genes          = dict()   # index the variants
    VariantAliases = dict()   # for possibly crosslinking databases
    Fasta          = dict()   # fasta sequences for reformatting variants

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
    parser.add_argument('-o', dest='output_file', type=str, required=False, help='output filename', default='druggability.out')
    parser.add_argument('-l', dest='log_file', type=str, required=False, help='logfile name', default='druggability.log')
    parser.add_argument('-d', '--debug', action='count', default=0)
    parser.add_argument('-nn', dest='normal_name', type=str, required=False, help='normal sample name', default='')
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-t', dest='variation_type', type=str, required=True, help='variation type:  maf | fusion')
    requiredNamed.add_argument('-f', dest='variant_file', type=str, required=True, help='variant filename')
    requiredNamed.add_argument('-tn', dest='tumor_name', type=str, required=True, help='tumor sample name')
    args = parser.parse_args()

    # Set up debug level
    if args.debug > 0:
        myglobal.DEBUG = True
    if args.debug > 1:
        myglobal.DEBUG_2 = True

    # Set up log file
    mylogfile = args.log_file
    if os.path.exists( mylogfile ):
        os.remove( mylogfile )
    logging.basicConfig( filename = mylogfile, level=logging.DEBUG)

    # Redirect output
    f = open( args.output_file, 'w')
    sys.stdout = f

    # Record run configuration
    logger = logging.getLogger('runconfig')
    logger.info('output file={}'.format( args.output_file ))
    logger.info('log file={}'.format( args.log_file ))
    logger.info('debug flag count={}'.format( args.debug ))
    logger.info('variant file={}'.format( args.variant_file ))
    logger.info('variation type={}'.format( args.variation_type ))
    logger.info('tumor sample name={}'.format( args.tumor_name ))
    logger.info('normal sample name={}'.format( args.normal_name ))
    logger.info('civic upstream version={}'.format( config.civic_params['upstream_version'] ))
    logger.info('civic liftover={}'.format( config.civic_params['ref_build_liftover'] ))

    # Validate input and run
    log_timestamp('start time')
    if args.variation_type not in ['maf', 'fusion']:
        abort_run('invalid input variation type')

    if args.variation_type == 'maf':
        if args.normal_name == '':
            abort_run('normal sample name is missing for maf input')

    if args.variation_type == 'fusion':
        if args.normal_name != '':
            logging.info('normal sample name will be ignored for fusion input')

    main( args )
    log_timestamp('end time')

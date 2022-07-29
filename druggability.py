#!/usr/bin/env python3

# R. Jay Mashl <rmashl@wustl.edu>

import os, sys, csv, re
import getopt
import argparse
import logging

import myglobal
from load_databases import *
from process_maf import *
from process_fusions import *
from utils import abort_run
import druggability_databases.config as config

# Update for module location
myglobal.DBPATH = os.path.join( os.path.dirname( os.path.abspath(__file__)), myglobal.DRUGDBDIR )


def main( args ):
    Matches        = dict()   # alteration matches by sample, separated in 'full' and 'partial' match lists
    Evidence       = dict()   # usage: key= drug, value = dict()
    Variants       = dict()   # variant records from alteration databases
    Genes          = dict()   # make list of variants for a given gene mentioned in alteration databases
    VariantAliases = dict()   # for possibly crosslinking databases
    Fasta          = dict()   # fasta sequences for reformatting variants
    Trials         = dict()   # clinical trials records
    Genes_altered  = dict()   # genes seen in input variant file
    Matches_trials = dict()   # trials matches by sample with breakdown by variant classes

    # load variant summary from CIViC
    load_civic( Variants, Genes, VariantAliases)

    # Load evidence from CIViC
    load_civic_evidence( Evidence, Variants )

    # Load fasta sequences
    load_fasta( Fasta )

    # Load clinical trials
    if len(args.annotate_trials):
        load_trials( Trials, args.annotate_trials )

    # call main processing
    if args.variation_type == 'maf':
        process_maf( args, Matches, Evidence, Variants, Genes, Fasta, Genes_altered, Trials, Matches_trials )

    if args.variation_type == 'fusion':
        process_fusions( args, Matches, Evidence, Variants, Genes, Genes_altered, Trials, Matches_trials )

    # Report results
    print_summary_for_all( args, Matches, Variants, Evidence, Matches_trials )


def list_trials():
    print('\nList of searchable clinical trials:\n')
    print( '\t'.join(['Keyword','Disease']) )
    print( '\t'.join(['-------','-------']))
    for k,v in config.trials_files.items():
        print( '\t'.join([ k, v['disease_name']]))
    print('\n')


if __name__ == '__main__':

    # Process command line
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', dest='output_file', type=str, required=False, help='output filename', default='druggability.out')
    parser.add_argument('-l', dest='log_file', type=str, required=False, help='logfile name', default='druggability.log')
    parser.add_argument('-d', '--debug', action='count', default=0)
    parser.add_argument('-nn', dest='normal_name', type=str, required=False, help='normal sample name', default='')
    parser.add_argument('-at', type=str, dest='annotate_trials', required=False, help='report clinical trials for this disease keyword', default='')
    parser.add_argument('-ato', dest='trials_auxiliary_output_file', type=str, required=False, help='clinical trials auxiliary output filename', default='trials.aux')
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-t', dest='variation_type', type=str, required=True, help='variation type:  maf | fusion')
    requiredNamed.add_argument('-f', dest='variant_file', type=str, required=True, help='variant filename')
    requiredNamed.add_argument('-tn', dest='tumor_name', type=str, required=True, help='sample or tumor sample name')
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
    logging.basicConfig( filename = mylogfile, level=0)

    # Validate input
    if args.variation_type not in ['maf', 'fusion']:
        abort_run('invalid input variation type')

    if args.variation_type == 'maf':
        if args.normal_name == '':
            abort_run('normal sample name is missing for maf input')

    if args.variation_type == 'fusion':
        if len(args.normal_name):
            logging.info('normal sample name will be ignored for fusion input')

    args.annotate_trials = args.annotate_trials.lower()
    if len(args.annotate_trials):
        if args.annotate_trials not in [ w.lower() for w in config.trials_files.keys()]:
            list_trials()
            abort_run('keyword ' + args.annotate_trials + ' does not have clinical trials annotations')

    # Redirect output
    f = open( args.output_file, 'w')
    sys.stdout = f

    # Record run configuration
    logger = logging.getLogger('runconfig')
    logger.setLevel(0)

    # Check version
    major,minor,micro = sys.version_info[0:3]
    logger.info('Python {}.{}.{} detected' . format(major,minor,micro))
    if (major != 3):
        abort_run('Python verion 3 is required')

    logger.info('command={}'.format( ' '.join( sys.argv )))
    logger.info('output file={}'.format( args.output_file ))
    logger.info('log file={}'.format( args.log_file ))
    logger.info('debug flag count={}'.format( args.debug ))
    logger.info('variant file={}'.format( args.variant_file ))
    logger.info('variation type={}'.format( args.variation_type ))
    if args.variation_type == 'maf':
        logger.info('tumor sample name={}'.format( args.tumor_name ))
        logger.info('normal sample name={}'.format( args.normal_name ))
    elif args.variation_type == 'fusion':
        logger.info('sample name={}'.format( args.tumor_name ))
    else:
        pass
    logger.info('civic upstream version={}'.format( config.civic_params['upstream_version'] ))
    logger.info('civic liftover={}'.format( config.civic_params['ref_build_liftover'] ))

    if len(args.annotate_trials):
        logger.info('clincal trials disease search={}'.format( args.annotate_trials ))
        logger.info('clincal trials parent accession={}'.format( config.trials_params['date_accessed'] ))

    # Run
    log_timestamp('start time')
    main( args )
    log_timestamp('end time')

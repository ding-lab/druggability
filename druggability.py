#!/usr/bin/env python3

# R. Jay Mashl <rmashl@wustl.edu>

import os, sys, csv, re
import getopt
import argparse
import logging
import pprint

import myglobal
from load_databases import *
from process_maf import *
from process_fusions import *
from utils import abort_run
import druggability_databases.config as config
from enums import *

# Update for module location
myglobal.DBPATH = os.path.join( os.path.dirname( os.path.abspath(__file__)), myglobal.DRUGDBDIR )


def log_sample_mentioned( SampleMentioned, key ):
    report_key = key
    if key == 'maf':
        report_key = 'maf/vcf'  # adjustment
    if SampleMentioned[ key ]:
        logger.info('Sample name was mentioned in the {} input file' . format(report_key))
    else:
        logger.info('WARNING: Sample name was NOT mentioned in the {} input file' . format(report_key))

def main( args ):
    Matches        = dict()   # alteration matches by sample, separated in 'full' and 'partial' match lists
    Evidence       = dict()   # usage: key= drug, value = dict()
    Variants       = dict()   # variant records from alteration databases
    Genes          = dict()   # make list of variants for a given gene mentioned in alteration databases
    VariantAliases = dict()   # for possibly crosslinking databases
    Fasta          = dict()   # fasta sequences for reformatting variants
    Trials         = dict()   # clinical trials records
    Genes_altered  = dict()   # genes seen in input variant file
    Matches_trials = dict()   # trials matches by sample with breakdown by variant classes; disqualified trials are a one-off key element
    Gene_sets      = dict()   # named gene sets/classes appearing trials
    GenesSeenInTrials = dict()
    SampleMentioned = dict()  # record whether sample was seen to QC against passed sample names

    # load variant summary from CIViC
    load_civic( Variants, Genes, VariantAliases)

    # Load evidence from CIViC
    load_civic_evidence( Evidence, Variants )

    # Load fasta sequences
    load_fasta( Fasta )

    # Load clinical trials
    if args.annotate_trials:
        load_gene_sets( Gene_sets )
        load_trials( Trials, args.annotate_trials, Gene_sets )
        if args.b_dump_trials_only:
            output_s = pprint.pformat( Trials, indent=1, width=100, sort_dicts=True)
            with open('trials.dump', 'w') as ff:
                ff.write( output_s )
            sys.exit(0)

    # Set up storage
    GenesSeenInTrials = dict()    #  Track genes present in the maf (by alteration type) that are relevant to trials in given call context
    for alt_type in [ MUTATION, INSERTION, DELETION, FUSION ]:
        GenesSeenInTrials[ alt_type ] = []
    GenesSeenInTrials['all_types'] = []    # one-off key to store merged list of all genes

    SampleMentioned = {'maf': False, 'fusion': False}

    if args.annotate_trials:
        # Use of a key for sample is legacy; use generic name to enable merging of modalities
        check_alloc_named( Matches_trials, SAMPLENAME, 'dict')
        for vc in VARIANT_CLASSES:
            check_alloc_named( Matches_trials[ SAMPLENAME ], vc, 'dict' )

    call_context = ''

    # call main processing with QC of sample names detected
    if args.variant_fusion_file:   # somatic fusions
        call_context = 'somatic'
        process_fusions( args, Matches, Evidence, Variants, Genes, Genes_altered, Trials, Matches_trials, SampleMentioned, call_context, GenesSeenInTrials )
        log_sample_mentioned( SampleMentioned, 'fusion')

    if args.varInputFile:       # somatic maf
        call_context = 'somatic'
        process_maf( args, Matches, Evidence, Variants, Genes, Fasta, Genes_altered, Trials, Matches_trials, SampleMentioned, call_context, GenesSeenInTrials )
        log_sample_mentioned( SampleMentioned, 'maf')  # or VCF

    # Process trials
    if args.annotate_trials:
        evaluate_gene_trials_hits( args, Trials, Genes_altered, SAMPLENAME, Matches_trials, call_context, GenesSeenInTrials )
        # Evaluate trial disqualifications
        Matches_trials[ SAMPLENAME ][ DISQUALIFYING ] = []
        if 'maf' in SampleMentioned.keys() or 'fusion' in SampleMentioned.keys():
            evaluate_trials_disqualifications( Matches_trials )


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
    parser.add_argument('-o', dest='output_file', type=str, required=False, help='alteration database matches', default='alterations.out')
    parser.add_argument('-l', dest='log_file', type=str, required=False, help='logfile name', default='druggability.log')
    parser.add_argument('-d', '--debug', action='count', default=0)
    parser.add_argument('-nn', dest='normal_name', type=str, required=False, help='MAF normal sample name', default='')
    parser.add_argument('-tn', dest='tumor_name', type=str, required=False, help='MAF tumor sample name', default='')
    parser.add_argument('-fn', dest='fusion_sample_name', type=str, required=False, help='fusion sample name', default='')
    parser.add_argument('-at', type=str, dest='annotate_trials', required=False, help='report clinical trials for this disease keyword', default='')
    parser.add_argument('-ato', dest='trials_auxiliary_output_file', type=str, required=False, help='clinical trials auxiliary output filename', default='trials.aux')
    parser.add_argument('--dump_trials_only', dest='b_dump_trials_only', action='store_true')
    parser.add_argument('--vcf', dest='variant_vcf_file', type=str, required=False, help='variant file in VCF format (requirements: tumor and normal names)', default='')
    parser.add_argument('--maf', dest='variant_maf_file', type=str, required=False, help='variant file in maf format (requirements: tumor and normal names)', default='')
    parser.add_argument('--basicmaf', dest='variant_basicmaf_file', type=str, required=False, help='variant file in basic maf format (requirements: tumor and normal names)', default='')
    parser.add_argument('--fusion', dest='variant_fusion_file', type=str, required=False, help='variant file for fusions (requirement: tumor name)', default='')
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
    if not len(sys.argv) > 1:
        parser.print_help()
    if not ( args.variant_vcf_file  or  args.variant_maf_file  or  args.variant_basicmaf_file  or  args.variant_fusion_file ):
        abort_run('please specify at least one file with alterations')

    varfileCount = 0
    bHasVarfile  = False
    varInputFile = ''
    for varfile in [ args.variant_maf_file, args.variant_basicmaf_file, args.variant_vcf_file ]:
        if varfile != '':
            varfileCount += 1
            args = namespace_append( args, 'varInputFile', varfile)   # save real filename
    if varfileCount > 1:
        abort_run('please specify at most one maf, basic maf, or VCF variant file')

    if bHasVarfile and not (args.normal_name and args.tumor_name):
        abort_run('tumor and normal sample names are required for maf or VCF input')

    if args.variant_fusion_file and not args.fusion_sample_name:
        abort_run('please specify fusion sample name')

    args.annotate_trials = args.annotate_trials.lower()
    if args.annotate_trials:
        if args.annotate_trials not in [ w.lower() for w in config.trials_files.keys()]:
            list_trials()
            abort_run('keyword ' + args.annotate_trials + ' does not have clinical trials annotations')


    # Record run configuration
    logger = logging.getLogger('runconfig')
    logger.setLevel(0)

    # Check version
    major,minor,micro = sys.version_info[0:3]
    logger.info('Python {}.{}.{} detected' . format(major,minor,micro))
    if (major < 3) or (major==3 and minor < 7) or (major==3 and minor==7 and micro < 6):
        abort_run('Python version 3.7.6 or higher is required')

    logger.info('command={}'.format( ' '.join( sys.argv )))
    logger.info('output file={}'.format( args.output_file ))
    logger.info('log file={}'.format( args.log_file ))
    logger.info('debug flag count={}'.format( args.debug ))

    if args.variant_vcf_file:
        logger.info('variant file={} (VCF format)'.format( args.variant_vcf_file ))
        logger.info('VCF tumor sample name={}'.format( args.tumor_name ))
        logger.info('VCF normal sample name={}'.format( args.normal_name ))
    if args.variant_maf_file:
        logger.info('variant file={} (maf format)'.format( args.variant_maf_file ))
        logger.info('MAF tumor sample name={}'.format( args.tumor_name ))
        logger.info('MAF normal sample name={}'.format( args.normal_name ))
    if args.variant_basicmaf_file:
        logger.info('variant file={} (basic maf format)'.format( args.variant_basicmaf_file ))
        logger.info('MAF tumor sample name={}'.format( args.tumor_name ))
        logger.info('MAF normal sample name={}'.format( args.normal_name ))
    if args.variant_fusion_file:
        logger.info('variant file={} (fusion format)'.format( args.variant_fusion_file ))
        logger.info('fusion sample name={}'.format( args.fusion_sample_name ))

    logger.info('civic upstream version={}'.format( config.civic_params['upstream_version'] ))
    logger.info('civic liftover={}'.format( config.civic_params['ref_build_liftover'] ))

    if args.annotate_trials:
        logger.info('clincal trials disease search={}'.format( args.annotate_trials ))
        logger.info('clincal trials access date={}'.format( config.trials_files[ args.annotate_trials ]['accessed'] ))
        logger.info('Dump trials only={}'.format( args.b_dump_trials_only ))

    # Run
    log_timestamp('start time')
    main( args )
    log_timestamp('end time')

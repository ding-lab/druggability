
# R. Jay Mashl <rmashl@wustl.edu>

import os, sys, csv, re
import config
from utils import *

DEBUG=config.DEBUG

def process_maf( inputFile, Evidence, Variants, Genes):

    Variant_sample_tracking = dict()   # record which samples have which variants

    tsv_file = open( inputFile )
    read_tsv = csv.reader(tsv_file, delimiter='\t')

    for row in read_tsv:
        fields = [ s.strip() for s in row ]

        # Header
        if re.search( r'^#', row[0]) or row[0] == "Hugo_Symbol":
            continue

        hugo         = fields[ 0]
        chrom        = fields[ 4]
        pos_start    = fields[ 5]
        pos_end      = fields[ 6]
        vartype      = fields[ 9]   # SNP, INS, DEL, ...
        sample       = fields[15]   # here, this is tumor sample
        aachange     = fields[36]
        csq          = fields[50]

        # ##################################################
        # initially look at only vars with AA change

        if not re.search( r'^p\.', aachange):
            continue

        aachange = aachange.replace('p.','')
        # ##################################################

        # Check whether this gene is mentioned in any database
        num_var_entries_for_gene_in_db = len(Genes[hugo]) if hugo in Genes.keys() else 0
        #if DEBUG:
        #    print( "Gene, num var entries for gene in db: " +  hugo + ', ' + str(num_var_entries_for_gene_in_db))
        if not num_var_entries_for_gene_in_db:
            continue

        # summarize alteration
        alteration_summary = '\t'.join([ hugo, chrom, pos_start, aachange, vartype ])


        # Merge lists of reported aliases
        dbSNP_RS                       = fields[ 13]  #rs#, list of rs#s, or novel
        existing_variation_list_string = fields[ 57]
        maf_variant_id                 = fields[111]
        known_variants = []
        if dbSNP_RS != '.':
            known_variants.extend( dbSNP_RS.split(','))
        if existing_variation_list_string != '.':
            known_variants.extend( existing_variation_list_string.split(','))
        if maf_variant_id != '.':
            known_variants.extend( maf_variant_id.split(','))
        known_variants = list(set( known_variants ))

        # Set up storage for tracking matches
        if sample not in Variant_sample_tracking.keys():
            Variant_sample_tracking[ sample ] = dict()
        Variant_sample_tracking[ sample ][alteration_summary] = dict( total_evidence_count=0, v_id_list=[] )



        # Identify and track matches
        num_v_id_matched     = 0
        num_v_id_not_matched = 0

        for v_id in Genes[hugo]:
            if is_exact_aa_match( aachange, v_id, Variants ) or is_pattern_aa_match( aachange, v_id, Variants):   # FIXME: false positives possible
                num_v_id_matched += 1
                #if DEBUG:
                #   print("Matched on (v_id, variant) = " + v_id + "," +  Variants[v_id]['variant'])
                Variant_sample_tracking[sample][alteration_summary]['v_id_list'].append( v_id )
                Variant_sample_tracking[sample][alteration_summary]['total_evidence_count'] += len(Variants[v_id]['evidence_list'])
            else:
                num_v_id_not_matched += 1

        #if num_v_id_matched > 0:
            #if DEBUG:
            #print('Summary: matched %s / not matched = %s of %s variants' % ( num_v_id_matched, num_v_id_not_matched, num_var_entries_for_gene_in_db))
            #print('')

    tsv_file.close()



    # Print summary
    for sample in Variant_sample_tracking.keys():
        for alteration in Variant_sample_tracking[sample].keys():
            this_alt = Variant_sample_tracking[sample][alteration]
            if this_alt['total_evidence_count']:
                print_sample_header( sample, alteration )
                print_output_header()
                for v_id in this_alt['v_id_list']:
                    for ev_id in Variants[v_id]['evidence_list']:
                        t = Evidence[ev_id]
                        print( *[ (v_id.split(':'))[0],   Variants[v_id]['variant'], t['disease'], t['oncogenicity'], t['mutation_effect'],   t['drugs_list_string'], t['evidence_type'], t['evidence_direction'], t['evidence_level'], t['clinical_significance'], format_citations(t['citations'])], sep = '\t')

                print('')
                print('')

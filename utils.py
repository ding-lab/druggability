
# R. Jay Mashl <rmashl@wustl.edu>

import os, sys, re
import druggability_databases.config as config
import logging
import datetime
from enums import *
import pandas as pd

logger = logging.getLogger(__name__)
logger.setLevel(0)

def abort_run( msg ):
    logger.error( msg )
    sys.exit(1)

def log_timestamp( label ):
    logger.info('{}: {}'.format(label, datetime.datetime.now() ))

# Convert list to string
def list2str( mylist ):
    return ','.join( filter( None, mylist ))

def is_exact_aachange_match( query, v_id, Variants ):
    if query == Variants[v_id]['variant']:
        return True
    else:
        return False

def has_genomic_match( query, v_id, Variants, fieldname ):
    if query == Variants[v_id][ fieldname ]:
        return True
    else:
        return False

def get_aachange_overlap_length( query, v_id, Variants ):
    a0, a1 = get_pos_range( query )
    b0, b1 = Variants[v_id]['prot_ref_start_pos'], Variants[v_id]['prot_ref_end_pos']

    if a0 <= b0:
        P0 = a0
        P1 = a1
        Q0 = b0
        Q1 = b1
    else:
        P0 = b0
        P1 = b1
        Q0 = a0
        Q1 = a1
    k = P1 - Q0 + 1
    return k if k > 0 else 0

def is_pattern_aa_match( aachange, v_id, Variants ):
    my_regex = r'^' + re.escape(Variants[v_id]['variant'])
    if re.search( my_regex, aachange, re.IGNORECASE):   # allows for missing alt aa but false hits are possible
        return True
    else:
        return False

# Decide sample display info
def get_sample_info_from_call( args, s ):
    if re.search( r'^MUT:', s ):
        sample_info = '||'.join([ args.tumor_name, args.normal_name ])
    elif re.search( r'^FUS:', s ):
        sample_info = args.fusion_sample_name
    else:
        abort_run('cannot determine sample description from call type')
    return sample_info

# Decide sample display info WILDTYPE, MUTATION, FUSION, INSERTION, DELETION
def get_sample_info_from_muttype( args, muttype ):
    if muttype in [ MUTATION, INSERTION, DELETION ]:
        sample_info = '||'.join([ args.tumor_name, args.normal_name ])
    elif muttype in [ FUSION ]:
        sample_info = args.fusion_sample_name
    elif muttype in [ WILDTYPE ]:
        if args.tumor_name and args.normal_name and not args.fusion_sample_name:
            sample_info = '||'.join([ args.tumor_name, args.normal_name ])
        elif args.fusion_sample_name and not (args.tumor_name or args.normal_name):
            sample_info = args.fusion_sample_name
        elif args.tumor_name and args.normal_name and args.fusion_sample_name:
            sample_info = '||'.join([ args.tumor_name, args.normal_name, args.fusion_sample_name ])
        else:
            abort_run('cannot construct sample description from sample names')
    else:
        abort_run('cannot determine sample description for trials output')
    return sample_info

def check_alloc_named( obj, key, s ):
    if key not in obj.keys():
        if s == 'list':
            obj[ key ] = []
        elif s == 'dict':
            obj[ key ] = dict()
        elif s == 'match_level':
            obj[ key ] = {'full': [], 'partial': [], 'wildtype': dict()}
        elif s == 'trial':
            obj[ key ] = dict()
            for i in VARIANT_CLASSES:
                obj[ key ][ i ] = dict()  # keys are the call contexts
        else:
            abort_run('Unexpected allocation type request')
    return

def format_citations( mylist ):
    result_list = []
    if len(mylist):
        for cit in mylist:
            if cit['source'].lower() == 'pubmed':
                result_list.append( '{source}:{myid}'.format( source='PMID', myid=cit['citation_id'] ) )
            else:
                result_list.append( '{source}:{myid}'.format( source=cit['source'], myid=cit['citation_id'] ) )
        return '; '.join( result_list )
    else:
        return '-'

def print_thin_line():
    print( '#' + '-' * 230 )

def print_thick_line():
    print( '#' * 70 )


header_merged_list    = [ 'Sample_info', 'Called',  'DB_Original', 'DB_Liftover', 'Match_Status', 'Criteria_Met',  'Source', 'Disease', 'Oncogenicity', 'Mutation_Effect', 'Treatment', 'Evidence_Type', 'Evidence_Direction',
                          'Evidence_Level', 'Clinical_Significance',  'Citation']
header_by_sample_list = [ 'Sample', 'Match_Index', 'Matched_Alteration', 'Match_Status', 'Criteria_Met',  'Source', 'Disease', 'Oncogenicity', 'Mutation_Effect', 'Treatment', 'Evidence_Type', 'Evidence_Direction',
                          'Evidence_Level', 'Clinical_Significance',  'Citation']
header_aux_list       = [ 'Sample_info', 'Disease', 'Variant_class', 'Call_context', 'Gene', 'Position_target', 'Trial_id', 'Intervention', 'Overall_status', 'Phase', 'Primary_completion_date']

def print_header( var_mode ):
    if var_mode in ['merged']:
        print( '\t'.join(header_merged_list) )
    elif var_mode == 'by_sample':   # currently unused
        print( '\t'.join(header_by_sample_list) )
    else:
        abort_run('unknown variant mode/filetype for printing')


def print_sample_header( sample, alteration ):
    print_thick_line()
    print( '### Sample:\t' + sample)
    print( '\t'.join([ '### Hugo_Symbol', 'Chromosome', 'Start_Position', 'AA_change', 'Variant_Type']))
    print( '### ' + alteration )
    print_thick_line()


def condense_altmatch_output( df, var_mode ):
    # remove copies
    df = df.drop_duplicates()

    # group citations
    if var_mode in ['merged']:
        header_orig    = header_merged_list.copy()
        header_groupby = header_merged_list.copy()
    else:
        abort_run('unknown variant mode/filetype in condense')
    header_groupby.remove('Citation')

    df = df.sort_values( by=df.columns.to_list() )
    df = df.groupby( header_groupby )['Citation'].apply(','.join).reset_index()
    df = df.reindex( columns = header_orig )
    df = df.sort_values( by=['Criteria_Met','Called'] )
    return df

def collapse_field( df, field ):
    header_aux_list_groupby = header_aux_list.copy()
    header_aux_list_groupby.remove( field )
    df = df.sort_values( by=[ field ] )
    df = df.groupby( header_aux_list_groupby )[field].apply(','.join).reset_index()
    df = df.reindex( columns = header_aux_list )
    return df

def condense_trials_output( df ):
    # remove copies
    df = df.drop_duplicates()

    # group genes, variant classes
    df = collapse_field( df, 'Gene' )
    df = collapse_field( df, 'Variant_class' )

    # combine genes in same family
    for idx, row in df.iterrows():
        genelist = row['Gene'].split(',')
        base_genes = dict()
        if len(genelist) > 1:
            for g in genelist:
                m = re.match(r'([A-Za-z]+)(\d+)?([A-Z]+)?([0-9]+)?', g)
                if m[2] is None:
                    check_alloc_named( base_genes, m[0], 'list' )
                    base_genes[ m[0] ].append( '' )
                elif m[3] is None:
                    check_alloc_named( base_genes, m[1], 'list' )
                    base_genes[ m[1] ].append( m[2] )
                elif m[4] is None:
                    check_alloc_named( base_genes, m[1]+m[2], 'list' )
                    base_genes[ m[1]+m[2] ].append( m[3] )    # resolve e.g. CA/CB later
                else:
                    check_alloc_named( base_genes, m[1]+m[2]+m[3], 'list' )
                    base_genes[ m[1]+m[2]+m[3] ].append( m[4] )
        else:
            check_alloc_named( base_genes, genelist[0], 'list' )
            base_genes[ genelist[0] ] = ''

        gene_families = []
        for gf in sorted(base_genes.keys()):
            gene_families.append( str(gf) + '/'.join(base_genes[gf]) )
        df.at[idx,'Gene'] = ','.join(gene_families)

    return df


def intersection( lst1, lst2 ):
    return list( set(lst1)  &  set(lst2) )

def uniquify( lst ):
    return list( set(lst) )

def list_append( lst, i ):
    if i not in lst:
        lst.append( i )
    return

# interpret protein positions
def parse_range( s, flag, infostr ):
    x = []                               # end points
    seg = s.split('_')
    if len(seg) > 2:
        abort_run('cannot identify range in ' + s)
    for si in seg:
        si = re.sub('[A-Y]', '', si)     # remove aa
        x.append( int(si) )
    if len(seg) == 1:
        x.append( int(x[0]) )        # make single position into interval
    return x

# get protein position range
def get_pos_range( v ):
    m = re.search(r'(.*?)(delins)(.*?)$', v)
    if m is not None:
        return parse_range( m[1], 'delins', m[3])

    m = re.search(r'(.*?)(del|ins|_splice)$', v)
    if m is not None:
        return parse_range( m[1], '', '')

    m = re.search(r'(.*?)(del|ins|dup)(.*?)$', v)
    if m is not None:
        return parse_range( m[1], '', '')

    m = re.search(r'(.*?)(fs)(.*?)$', v)
    if m is not None:
        return parse_range( m[1], 'fs', m[3])

    m = re.search(r'(.*?)(ext\*|ext)(.*?)$', v)
    if m is not None:
        return parse_range( m[1], 'ext', m[3])

    m = re.search(r'([A-Y])(\d+)([A-Y]|\*)?$', v)
    if m is not None:
        return parse_range( m[2], '', '')

    if v == 'M1?':
        return parse_range( 'M1', '', '')

    abort_run('unable to determine position range for ' + v)


# calculate genomic coordinate range
# ...a stop coordinate is ld be generated
def calculate_gdna_coords( variant_set, liftover_status ):
    if liftover_status == 'use_liftover':
        pos0 = variant_set['start_liftover']
        pos1 = variant_set['stop_liftover']
    else:
        pos0 = variant_set['pos0']
        pos1 = variant_set['pos1']

    if not len(pos1):
        pos1 = pos0     # create stop coordinate where not provided; needed for checking interval overlap
    return '{chrom}:g.{pos0}_{pos1}'.format( chrom=variant_set['chrom'], pos0=pos0, pos1=pos1 )


# Genomic overlap with liftover coords
def has_genomic_overlap( query_dict, v_id, Variants ):
    q_chrom = query_dict['chrom']
    if q_chrom != Variants[v_id]['chrom_liftover']:
        return int(0)
    a0, a1 = int(query_dict['pos0']), int(query_dict['pos1'])
    b0, b1 = int(Variants[v_id]['start_liftover']), int(Variants[v_id]['stop_liftover'])

    if a0 <= b0:
        P0 = a0
        P1 = a1
        Q0 = b0
        Q1 = b1
    else:
        P0 = b0
        P1 = b1
        Q0 = a0
        Q1 = a1
    k = P1 - Q0 + 1
    return k if k > 0 else 0


def calculate_gdna_change( variant_set, liftover_status ):
    gene      = variant_set['gene']
    ref       = variant_set['ref']
    alt       = variant_set['alt']

    if liftover_status == 'use_liftover':
        chrom     = variant_set['chrom_liftover']
        start_pos = variant_set['start_liftover']
        end_pos   = variant_set['stop_liftover']
    else:
        chrom     = variant_set['chrom']
        start_pos = variant_set['pos0']
        end_pos   = variant_set['pos1']


    # Determine mutation type (exclude fusions to avoid adding 'delins' attribute
    if variant_set['main_variant_class'] == FUSION:
        return '{chrom}:g.{start}_{stop}'.format( chrom=chrom, start=start_pos, stop=end_pos )

    if (not len(ref)) and (    len(alt)):
        if (int(end_pos) - int(start_pos)) != 1:
            logger.info('(TODO) positions may describe duplication rather than insertion ({gene} at g.{pos})'.format( gene=gene, pos=start_pos ))
            return ''
        return '{chrom}:g.{start}_{stop}ins{alt}'.format( chrom=chrom, start=start_pos, stop=end_pos, alt=alt )

    if (    len(ref)) and (not len(alt)):
        return '{chrom}:g.{start}_{stop}del' .format( chrom=chrom, start=start_pos, stop=end_pos )

    if len(ref) == 1 and len(alt) == 1:   # substitution
        return '{chrom}:g.{start}{ref}>{alt}'.format( chrom=chrom, start=start_pos, ref=ref, alt=alt )
    else:
        return '{chrom}:g.{start}_{stop}delins{alt}'.format( chrom=chrom, start=start_pos, stop=end_pos, alt=alt )


def print_summary_by_sample( Variant_tracking, Variants, Evidence ):
    for sample in Variant_tracking.keys():
            for alteration in Variant_tracking[sample].keys():
                this_alt = Variant_tracking[sample][alteration]
                if this_alt['total_evidence_count']:
                    print_sample_header( sample, alteration )
                    print_header('by_sample')
                    for v_id in this_alt['v_id_list']:
                        for ev_id in Variants[v_id]['evidence_list']:
                            t = Evidence[ev_id]
                            print( *[ v_id.split(':')[0],   Variants[v_id]['variant'], t['disease'], t['oncogenicity'], t['mutation_effect'],   t['drugs_list_string'], t['evidence_type'], t['evidence_direction'], config.evidence_level_anno[t['evidence_level']], t['clinical_significance'], format_citations(t['citations'])], sep = '\t')

                    print('')
                    print('')


# given a list as a string, return an array of the items
def clean_split( s ):
    return [i.strip() for i in re.split(',|;', s)]  # semicolon can appear in maf input

# map mutation type
def map_mut( s ):
    the_map = {'mutation': MUTATION, 'insertion': INSERTION, 'deletion': DELETION, 'fusion': FUSION, 'none': WILDTYPE}
    return the_map[s]

def map_mut_reverse( s ):
    the_map = {MUTATION: 'mutation', INSERTION: 'insertion', DELETION: 'deletion', FUSION: 'fusion', WILDTYPE: 'wildtype'}
    return the_map[s]

def map_maf_reverse( s ):
    the_map = {WASHU_MAF: 'WASHU_MAF', UNION_MAF: 'UNION_MAF', BASIC_MAF: 'BASIC_MAF'}
    return the_map[s]

def list_disqualified_trials( disquals ):
    mylist = []
    for i in disquals:
        mylist.append( disquals['trial_id'] )
    return uniquify( mylist )

def print_summary_for_all( args, Matches, Variants, Evidence, Matches_trials ):
    for s in Matches:                  # legacy loop over samples
        altmatch_output_lines = []     # output for raw matches to "alt"erations

        # Check whether there are results
        num_results = 0
        for matchtype in ['full', 'partial']:
            num_results += len(Matches[s][matchtype])

        if num_results == 0:
            orig_stdout = sys.stdout
            with open( args.output_file, 'w') as f:
                sys.stdout = f
                print_header( 'merged' if (args.variant_maf_file or args.variant_basicmaf_file or args.variant_fusion_file) else '' )
                sys.stdout = orig_stdout
                logger.info('No matches to alteration db found!')

        else:
            match_idx = 0
            for matchtype in ['full', 'partial']:
                if len(Matches[s][matchtype]):
                    for dic in Matches[s][matchtype]:
                        v_id   = dic['v_id']
                        reason = dic['reason']
                        called = dic['called']
                        for ev_id in Variants[v_id]['evidence_list']:
                            t = Evidence[ev_id]
                            match_idx += 1
                            db_orig_str     = '{v_id}|{variant}|{gchange}|{refbuild}'.format( v_id=v_id, variant=Variants[v_id]['variant'], gchange=Variants[v_id]['gdnachange'], refbuild=Variants[v_id]['ref_build'] )
                            db_liftover_str = '{v_id}|{variant}|{gchange}|{refbuild}'.format( v_id=v_id, variant=Variants[v_id]['variant'], gchange=Variants[v_id]['gdnachange_liftover'], refbuild=Variants[v_id]['ref_build_liftover'] )
                            sample_info     = get_sample_info_from_call( args, called )
                            altmatch_output_lines.append( [ sample_info, called, db_orig_str, db_liftover_str, matchtype, reason,   v_id.split(':')[0],  t['disease'], t['oncogenicity'], t['mutation_effect'],   t['drugs_list_string'], t['evidence_type'], t['evidence_direction'], config.evidence_level_anno[t['evidence_level']], t['clinical_significance'], format_citations(t['citations'])] )

            # use pandas to prepare output
            df = pd.DataFrame( altmatch_output_lines, columns = header_merged_list )
            df = condense_altmatch_output( df, 'merged' )
            df.to_csv( args.output_file, sep = '\t', header=True, index=False)

        # output the matches to trials
        if len(args.annotate_trials):
            disqualified_trials_list = list_disqualified_trials( Matches_trials[s][DISQUALIFYING] )
            logger.info('Number of clinical trial disqualifying alterations: {}'.format( len(Matches_trials[s][ DISQUALIFYING ])))
            logger.info('Trials disregarded due to disqualifying alterations: {}'.format( 'none' if not disqualified_trials_list else ','.join( disqualified_trials_list ) ))

            aux_output_lines = []
            for vt in VARIANT_CLASSES:   # reminder: wt, ins, del, mut, fusion
                if len( Matches_trials[s][ vt ] ):
                    for gene in Matches_trials[s][ vt ]:
                        for ct_info in Matches_trials[s][ vt ][ gene ]:
                            if ct_info['trial_id'] not in disqualified_trials_list:
                                sample_info  = get_sample_info_from_muttype( args, vt )
                                aux_output_lines.append([ sample_info, args.annotate_trials, map_mut_reverse( vt ), ct_info['call_context'], gene, ct_info['position_target'], ct_info['trial_id'], ct_info['intervention'], ct_info['overall_status'], ct_info['phase'], ct_info['completion_date']])

            # use pandas to prepare output
            df = pd.DataFrame( aux_output_lines, columns = header_aux_list )
            df = condense_trials_output( df )
            df.to_csv( args.trials_auxiliary_output_file, sep = '\t', header=True, index=False)

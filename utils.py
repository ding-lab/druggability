
# R. Jay Mashl <rmashl@wustl.edu>

import re

# Convert list to string
def list2str( mylist ):
    return ','.join( filter( None, mylist ))

def is_exact_aa_match( aachange, v_id, Variants ):
    if aachange == Variants[v_id]['variant']:
        return True
    else:
        return False

def is_pattern_aa_match( aachange, v_id, Variants ):
    my_regex = r'^' + re.escape(Variants[v_id]['variant'])
    if re.search( my_regex, aachange, re.IGNORECASE):   # allows for missing alt aa but false hits are possible
        return True
    else:
        return False

def format_citations( mylist ):
    result_list = []
    if len(mylist):
        for cit in mylist:
            if cit['source'].lower() == 'pubmed':
                result_list.append( 'PMID' + ':' + cit['citation_id'])
            else:
                result_list.append( cit['source'] + ':' + cit['citation_id'])
        return '; '.join( result_list )
    else:
        return '-'

def print_thin_line():
    print( '#' + '-' * 230 )

def print_thick_line():
    print( '#' * 70 )

def print_output_header():
    print( '\t'.join([ '#Sample', 'Match_Index', 'Matched_Alteration', 'Match_Status', 'Unchecked_Criteria',  'Source', 'Disease', 'Oncogenicity', 'Mutation_Effect', 'Treatment', 'Evidence_Type', 'Evidence_Direction',
           'Evidence_Level', 'Clinical_Significance',  'Citation']) )
    print_thin_line()

def print_sample_header( sample, alteration ):
    print_thick_line()
    print( '### Sample:\t' + sample)
    print( '\t'.join([ '### Hugo_Symbol', 'Chromosome', 'Start_Position', 'AA_change', 'Variant_Type']))
    print( '### ' + alteration )
    print_thick_line()

def intersection( lst1, lst2 ):
    return list( set(lst1)  &  set(lst2) )

def uniquify( lst ):
    return list( set(lst) )

def list_append( lst, i ):
    if i not in lst:
        lst.append( i )
    return

def print_summary_by_sample( Variant_tracking, Variants, Evidence ):
    for sample in Variant_tracking.keys():
            for alteration in Variant_tracking[sample].keys():
                this_alt = Variant_tracking[sample][alteration]
                if this_alt['total_evidence_count']:
                    print_sample_header( sample, alteration )
                    print_output_header()
                    for v_id in this_alt['v_id_list']:
                        for ev_id in Variants[v_id]['evidence_list']:
                            t = Evidence[ev_id]
                            print( *[ v_id.split(':')[0],   Variants[v_id]['variant'], t['disease'], t['oncogenicity'], t['mutation_effect'],   t['drugs_list_string'], t['evidence_type'], t['evidence_direction'], t['evidence_level'], t['clinical_significance'], format_citations(t['citations'])], sep = '\t')

                    print('')
                    print('')


def print_summary_for_all( Matches, Variants, Evidence ):
    bPrintHeader = True
    for s in Matches:
        match_idx = 0
        for matchtype in ['full', 'partial']:
            if len(Matches[s][matchtype]):
                for dic in Matches[s][matchtype]:
                    v_id   = dic['v_id']
                    reason = dic['reason']
                    for ev_id in Variants[v_id]['evidence_list']:
                        t = Evidence[ev_id]
                        match_idx += 1
                        if bPrintHeader:
                            print_output_header()
                            bPrintHeader = False
                        print( *[ s, match_idx,  Variants[v_id]['variant'], matchtype, reason,   v_id.split(':')[0],  t['disease'], t['oncogenicity'], t['mutation_effect'],   t['drugs_list_string'], t['evidence_type'], t['evidence_direction'], t['evidence_level'], t['clinical_significance'], format_citations(t['citations'])], sep = '\t')

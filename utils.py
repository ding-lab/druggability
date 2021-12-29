
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
    result = ""
    for cit in mylist:
        if result != "":
            result = result + '; '
        if cit['source'].lower() == 'pubmed':
            result = result + 'PMID' + ':' + cit['citation_id']
        else:
            result = result + cit['source'] + ':' + cit['citation_id']
    return result

def print_thin_line():
    print( '###' + '-' * 67 )

def print_thick_line():
    print( '#' * 70 )

def print_output_header():
    print( '\t'.join([ '### Source', 'Evidence_Variant', 'Disease', 'Oncogenicity', 'Mutation_Effect', 'Treatments', 'Evidence_Type', 'Evidence_Direction',
           'Evidence_Level', 'Clinical_Significance',  'Citations']) )
    print_thin_line()

def print_sample_header( sample, alteration ):
    print_thick_line()
    print( '### Sample:\t' + sample)
    print( '\t'.join([ '### Hugo_Symbol', 'Chromosome', 'Start_Position', 'AA_change', 'Variant_Type']))
    print( '### ' + alteration )
    print_thick_line()

def intersection( lst1, lst2 ):
    return list( set(lst1)  &  set(lst2) )

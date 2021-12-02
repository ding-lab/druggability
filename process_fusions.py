
# R. Jay Mashl <rmashl@wustl.edu>

import os, sys, csv, re
from utils import *

def process_fusions( inputFile, Evidence, Variants, Genes):
    tsv_file = open( inputFile )
    read_tsv = csv.reader(tsv_file, delimiter='\t')

    col = []    # field names
    for row in read_tsv:
        fields = [ s.strip() for s in row ]

        # Header
        if re.search( r'^#', row[0]) or row[0] == "FusionName":
            continue

        FusionName, LeftBreakpoint, RightBreakpoint, Sample, JunctionReadCount, SpanningFragCount, FFPM, PROT_FUSION_TYPE, GTEx, CallerN = fields


    ## WORK IN PROGRESS

    tsv_file.close()


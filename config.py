
DEBUG=True
#DEBUG=False

#DEBUG_2=True
DEBUG_2=False

druggability_databases_gitpath = '/Users/rmashl/git/druggability_databases'

# Read preprocessed DrugBank
drugbank_files = {
    'file':  ''
}

# Read in CIViC
civic_files = {
    'assertions':    druggability_databases_gitpath + '/CIViC/' + '01-Oct-2021-AssertionSummaries.tsv',
    'genes':         druggability_databases_gitpath + '/CIViC/' + '01-Oct-2021-GeneSummaries.tsv',
    'variants':      druggability_databases_gitpath + '/CIViC/' + '01-Oct-2021-VariantSummaries.tsv',
    'evidence':      druggability_databases_gitpath + '/CIViC/' + '01-Oct-2021-ClinicalEvidenceSummaries.tsv',
    'variantgroups': druggability_databases_gitpath + '/CIViC/' + '01-Oct-2021-VariantGroupSummaries.tsv',
    'variants_preprocessed':      druggability_databases_gitpath + '/CIViC/' + '01-Oct-2021-VariantSummaries.preprocessed.tsv',
}

# Read in current oncokb
oncokb_files = {
    'variants':      druggability_databases_gitpath + '/OncoKB/' + 'oncokb.annotated.tsv',
    'therapeutics':  druggability_databases_gitpath + '/OncoKB/' + 'oncokb.therapeutic.tsv',
}

# Read in gzipped fasta protein file
uniprot_files = {
    'fasta':         druggability_databases_gitpath + '/Uniprot/' + 'uniprot_sprot.fasta.human.tsv.gz',
}

evidence_level_anno = {
    # CIViC
    'A': 'approved', # renamed from /validated/
    'B': 'clinical',
    'C': 'case study',
    'D': 'preclinical',
    'E': 'inferential',

    # OncoKB
    '1': 'approved', # is: FDA-recognized biomarker predictive of response to an FDA-approved drug in this indication
    '2': 'approved', # is: standard care biomarker recommended by the NCCN or other professional guidelines predictive of response to an FDA-approved drug in this indication
    '3A': 'clinical', # is: compelling clinical evidence supports the biomarker as being predictive of repsonse to a drug in this indication
    '3B': 'clinical', # is: standard care or investigational biomarker predictive of response to an FDA-approved or investigational drug in another indiation
    '4': 'preclinical', # is: compelling biological evidence supports the biomarker as being predictive of response to a drug
    'Dx1': 'required for diagnosis', # is: FDA and/or professional guideline-recognized biomarker required for diagnosis in this indication
    'Dx2': 'supports diagnosis', # is: FDA and/or professional guildeline-recognized biomarker tha supports diagnosis in this indication
    'Dx3': 'assists diagnosis', # is: biomarker that may assist disease diagnosis in this indication based on clinical evidence
    'R1': 'resistance, standard care', # is: standard care biomarker predictive of resistance to an FDA-approved drug in this indication
    'R2': 'resistance, investigational', # is: compelling clinical evidence supports the biomarker as being predictive of resistance to a drug

    # Miscellaneous
    '-': '-',
}

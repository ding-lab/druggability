
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

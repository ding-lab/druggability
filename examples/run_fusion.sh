#!/bin/bash

# main script
DRUGGABILITY=../druggability.py

# Provide the location of a test maf; here, this is a cancer cohort
cohort_fusions=../Test_files/CPTAC_fusions_v0.1.csv.tsv

# Get samples
tail -n +2 $cohort_fusions | cut -f6 | sort -u > fusion_samples.tsv

OUTDIR=run_somatic_cohort_fusion_output
mkdir -p $OUTDIR

for tumor_sample in $(cat fusion_samples.tsv) ; do
    fileout_stem=$(basename $cohort_fusions).${tumor_sample}
    echo $fileout_stem
    $DRUGGABILITY  -t fusion -f $cohort_fusions -tn ${tumor_sample}  -l $OUTDIR/$fileout_stem.log  -o $OUTDIR/$fileout_stem.out
done

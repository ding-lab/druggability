#!/bin/bash

# main script. Druggability databases are a submodule called druggability_databases
DRUGGABILITY=~/git/druggability/druggability.py

# Provide the location of a test maf; here, this is a cancer cohort
cohort_maf=../Test_files/LUAD.Somatic.050919.mnp.annot.maf


# Get samples
cut -f16,17 $cohort_maf | sort -u | grep -v -i barcode > case_pairs.tsv

cat case_pairs.tsv | while read line ; do
    read tumor_sample normal_sample <<< $(echo "$line" | cut -d$'\t' -f 1,2)
    OUTDIR=run_somatic_cohort_output
    mkdir -p $OUTDIR
    fileout_stem=$(basename $cohort_maf).${tumor_sample}_${normal_sample}
    echo $fileout_stem
    $DRUGGABILITY  -t maf  -f  $cohort_maf  -nn ${normal_sample} -tn ${tumor_sample}  -l $OUTDIR/$fileout_stem.log  -o $OUTDIR/$fileout_stem.out
done

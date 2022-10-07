#!/bin/bash

# main script. Druggability databases are a submodule called druggability_databases
DRUGGABILITY=../druggability.py

# Provide the location of a test maf; here, this is a cancer cohort
cohort_maf=../Test_files/chol_paper.maf

# case pairs file
case_pairs_file=cases.3.tsv

# Get samples
# The pseudo maf was constructed with tumor = normal = case)
tail -n +2 $cohort_maf |cut -f10,11 | sort -u > $case_pairs_file

OUTDIR=run_somatic_chol_paper_for_trials_testing
mkdir -p $OUTDIR

cat $case_pairs_file | while read line ; do
    read tumor_sample normal_sample <<< $(echo "$line" | cut -d$'\t' -f 1,2)
    fileout_stem=chol_paper.${tumor_sample}_${normal_sample}
    echo $fileout_stem
    $DRUGGABILITY  -t basicmaf -f  $cohort_maf -nn ${normal_sample} -tn ${tumor_sample}  -l $OUTDIR/$fileout_stem.log  -o $OUTDIR/$fileout_stem.out  -at chol -ato $OUTDIR/$fileout_stem.aux -d
done

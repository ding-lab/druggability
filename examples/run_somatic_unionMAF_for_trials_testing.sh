#!/bin/bash

# main script. Druggability databases are a submodule called druggability_databases
DRUGGABILITY=~/git/druggability/druggability.py

# Provide the location of a test maf; here, this is a cancer cohort
cohort_maf=../Test_files/PanCan_Union_Maf_Broad_WashU_v1.1.maf

# case pairs file
case_pairs_file=case_pairs.1.tsv

# Get samples
tail -n +2 $cohort_maf | cut -f18,19 | sort -u > $case_pairs_file

cat $case_pairs_file | while read line ; do
    read tumor_sample normal_sample <<< $(echo "$line" | cut -d$'\t' -f 1,2)
    OUTDIR=run_somatic_unionMAF_output_for_trials_testing
    mkdir -p $OUTDIR
    fileout_stem=$(basename $cohort_maf).${tumor_sample}_${normal_sample}
    echo $fileout_stem
    $DRUGGABILITY  -t maf  -f  $cohort_maf  -nn ${normal_sample} -tn ${tumor_sample}  -l $OUTDIR/$fileout_stem.log  -o $OUTDIR/$fileout_stem.out   -at chol -ato $OUTDIR/$fileout_stem.aux -d
done

#!/bin/bash

# 2023-04-12: perform joint search

# main script. Druggability databases are a submodule called druggability_databases
DRUGGABILITY=../druggability.py

# Provide the location of test files (not distributed here)
cohort_maf=../Test_files/mc3.v0.2.8.PUBLIC.code.filtered.Missense.maf.CHOL
cohort_fusions=../Test_files/TCGA_MC3.allfusion.compile.txt.CHOL.reformatted.tsv

# Case/sample info: columns are Subject, MAF_tumor_name,  MAF_normal_name,  fusion_sample_name
case_pairs_file=../Test_files/MC3_sample_map.tsv

OUTDIR=run_somatic_TCGA_MC3_chol_for_trials.merged_code.mut_only
mkdir -p $OUTDIR
tail -n +2  $case_pairs_file | while read line ; do
    read subjID tumor_sample normal_sample fusion_sample <<< $(echo "$line" | cut -d$'\t' -f 1-4)
    fileout_stem=TCGA_MC3_chol.${subjID}
    echo $fileout_stem
    $DRUGGABILITY  --basicmaf  $cohort_maf   -nn ${normal_sample} -tn ${tumor_sample}  -l $OUTDIR/$fileout_stem.log  -o $OUTDIR/$fileout_stem.out  -at chol -ato $OUTDIR/$fileout_stem.aux -d
done

OUTDIR=run_somatic_TCGA_MC3_chol_for_trials.merged_code.fus_only
mkdir -p $OUTDIR
tail -n +2  $case_pairs_file | while read line ; do
    read subjID tumor_sample normal_sample fusion_sample <<< $(echo "$line" | cut -d$'\t' -f 1-4)
    fileout_stem=TCGA_MC3_chol.${subjID}
    echo $fileout_stem
    $DRUGGABILITY   --fusion $cohort_fusions   -fn ${fusion_sample} -l $OUTDIR/$fileout_stem.log  -o $OUTDIR/$fileout_stem.out  -at chol -ato $OUTDIR/$fileout_stem.aux -d
done

OUTDIR=run_somatic_TCGA_MC3_chol_for_trials.merged_code.both
mkdir -p $OUTDIR
tail -n +2  $case_pairs_file | while read line ; do
    read subjID tumor_sample normal_sample fusion_sample <<< $(echo "$line" | cut -d$'\t' -f 1-4)
    fileout_stem=TCGA_MC3_chol.${subjID}
    echo $fileout_stem
    $DRUGGABILITY  --basicmaf  $cohort_maf  --fusion $cohort_fusions  -nn ${normal_sample} -tn ${tumor_sample} -fn ${fusion_sample} -l $OUTDIR/$fileout_stem.log  -o $OUTDIR/$fileout_stem.out  -at chol -ato $OUTDIR/$fileout_stem.aux -d
done

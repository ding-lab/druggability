# druggability

## Requirements
- Software environment with python v3.7.6 or higher and pandas v1.0.5 or higher
- Alteration databases
- Variant file (single, or multi-sample, MAF or fusion results file; VCF not yet supported) aligned to GRCh38
- (Optional) Clinical trials metadata file(s)

## Installation
#### Local installation
1. A software environment such as `conda` or **[miniconda](https://docs.conda.io/en/latest/miniconda.html)** is highly recommended.
2. Get the repositories:

   `git clone --recurse-submodules https://github.com/ding-lab/druggability.git`

   Running this command creates the main `druggability` directory and the subdirectory `druggability_databases` with the database files. (It is possible to clone the repos separately; just be sure to update `DRUGDBDIR` in `myglobal.py`.)

3. If relevant, activate the conda/miniconda environment and install additional packages: `pip install -r requirements.txt`.

#### Docker environment
1. Obtain the `Dockerfile` by direct download or from a "git clone" of the `druggability` repo and then build the image, e.g.,

   `docker build -t  druggability .`


## Running
```
$ ./druggability.py -h
usage: druggability.py [-h] [-o OUTPUT_FILE] [-l LOG_FILE] [-d] [-nn NORMAL_NAME]
                       [-tn TUMOR_NAME] [-fn FUSION_SAMPLE_NAME] [-at ANNOTATE_TRIALS]
                       [-ato TRIALS_AUXILIARY_OUTPUT_FILE] [--dump_trials_only]
                       [--maf VARIANT_MAF_FILE] [--basicmaf VARIANT_BASICMAF_FILE]
                       [--fusion VARIANT_FUSION_FILE]

optional arguments:
  -h, --help                         show this help message and exit
  -o OUTPUT_FILE                     alteration database matches
  -l LOG_FILE                        logfile name
  -nn NORMAL_NAME                    MAF normal sample name
  -tn TUMOR_NAME                     MAF tumor sample name
  -fn FUSION_SAMPLE_NAME             fusion sample name
  -at ANNOTATE_TRIALS                report clinical trials for this disease keyword
  -ato TRIALS_AUXILIARY_OUTPUT_FILE  clinical trials auxiliary output filename
  --maf VARIANT_MAF_FILE             variant file in maf format (requirements: tumor and normal names)
  --basicmaf VARIANT_BASICMAF_FILE   variant file in basic maf format (requirements: tumor and normal names)
  --fusion VARIANT_FUSION_FILE       variant file for fusions (requirement: tumor name)
  --dump_trials_only
  -d, --debug
```
Notes:
- Releases v1.2+:
   - Maf and fusion files can now be specified at the same time. The sample(s) responsible for the database matches are shown in the output tables.
   - Each input file type has been given its own command-line flag.
   - The sample names passed to the tool (`-nn`, `-tn`, `-fn`) need to match those appearing in the inputs. This feature can be helpful for comparing results across aliquots.
- For fusion inputs, the "fusion sample name" can be any sample, normal or tumor.
- See the `examples` directory for driver scripts to use as templates, particularly for analyzing a cohort of sample sets.
- The variant call files that were used in development are not able to be distributed externally.

### MAF examples:
1. **MAFs** in ding-lab or PanCan union maf formats are automatically detected:
    ```
    $ ./druggability.py  --maf  LUAD.Somatic.050919.mnp.annot.maf -nn C3N-00169_N -tn C3N-00169_T  -l  druggability.log  -o alterations.out
    ```
    ```
    $ ./druggability.py  --maf  LUAD.Somatic.050919.mnp.annot.maf -nn C3N-00169_N -tn C3N-00169_T  -l druggability.log  -o alterations.out  -at chol  -ato  trials.aux
    ```

2. **Ad hoc MAFs**. Variant call files, such as those redistributed with publications, may have a different format than those listed above. In that case, `druggability` can read a tab-delimited file whose first 13 fields are the following:
    ```
    Hugo_Symbol	NCBI_Build	Chromosome	Start_Position	End_Position	Variant_Classification	Variant_Type	Reference_Allele	Tumor_Seq_Allele1	Tumor_Sample_Barcode	Matched_Norm_Sample_Barcode	HGVSc	HGVSp_Short
    ```

    and analyzed with the `basicmaf` option:

    ```
    $ ./druggability.py  --basicmaf  cohort.maf -nn 111 -tn 111  -l druggability.log  -o alterations.out  -at chol  -ato  trials.aux
    ```
    Notes:
	- Coordinates are 1-based
	- `Variant_Type` should be SNP, DNP, MNP, INS, or DEL

### Fusion example:
- Input file must be in ding-lab format.
    ```
    $ ./druggability.py  --fusion  CPTAC_fusions_v0.1.csv.tsv -fn C3L-00165_T  -l druggability.log  -o alterations.out
    ```

### Joint MAF/Fusion example:

    ```
    $ ./druggability.py --basicmaf mc3.v0.2.8.PUBLIC.code.filtered.Missense.maf.CHOL --fusion TCGA_MC3.allfusion.compile.txt.CHOL.reformatted.tsv -nn TCGA-3X-AAVA-10A-01D-A41A-09 -tn TCGA-3X-AAVA-01A-11D-A417-09 -fn TCGA-3X-AAVA-01A-11R-A41I-07 -l TCGA_MC3_chol.TCGA-3X-AAVA.log -o TCGA_MC3_chol.TCGA-3X-AAVA.out -at chol -ato TCGA_MC3_chol.TCGA-3X-AAVA.aux
    ```

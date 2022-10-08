# druggability

## Requirements
- Software environment with python v3.7.6 or higher and pandas v1.0.5 or higher
- Alteration databases
- Variant file (single, or multi-sample, MAF or fusion results file; VCF not yet supported) aligned to GRCh38
- (Optional) Clinical trials metadata file(s)

## Installation
#### Local installation
1. A software environment such as `conda` or **[miniconda](https://docs.conda.io/en/latest/miniconda.html)** is highly recommended if it is not otherwise provided for.
2. Get the repositories. Using SSH URL method,
	- `git clone git@github.com:ding-lab/druggability.git `
	- `cd druggability`
	- `git clone git@github.com:ding-lab/druggability_databases.git`
	
	In theory, a single `git clone --recurse-submodules ...` command should work once the repos become public.
3. If appropriate, activate the conda/miniconda environment and install additional packages: `pip install -r requirements.txt`.
4. Obtain our clinical trials database. As these files are being updated and tested continually, they are not being distributed externally at this time; please contact the author for availability status.  In the [druggability_databases](https://github.com/ding-lab/druggability_databases/) repository, edit `config.py` with the correct paths.

#### Docker environment



## Running
```
$ ./druggability.py -h
usage: druggability.py [-h] [-o OUTPUT_FILE] [-l LOG_FILE] [-d]
                       [-nn NORMAL_NAME] [-at ANNOTATE_TRIALS]
                       [-ato TRIALS_AUXILIARY_OUTPUT_FILE] -t VARIATION_TYPE
                       -f VARIANT_FILE -tn TUMOR_NAME

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT_FILE        alteration database matches
  -l LOG_FILE           logfile name
  -d, --debug
  -nn NORMAL_NAME       normal sample name
  -at ANNOTATE_TRIALS   report clinical trials for this disease keyword
  -ato TRIALS_AUXILIARY_OUTPUT_FILE
                        clinical trials auxiliary output filename

required arguments:
  -t VARIATION_TYPE     variation type: maf | fusion | basicmaf
  -f VARIANT_FILE       variant filename
  -tn TUMOR_NAME        sample or tumor sample name
```
Notes:
- The code is developed with cohort input files in mind and so it's best run as part of a loop. See the `examples` directory for driver scripts to use as templates.
- The variant call files that were used in development are not able to be distributed externally.
- For fusion inputs, the "tumor name" is thought of as any sample, normal or tumor.

### MAF examples:
1. **MAFs** in ding-lab or PanCan union maf formats are automatically detected:
    ```
    $ ./druggability.py  -t maf  -f LUAD.Somatic.050919.mnp.annot.maf -nn C3N-00169_N -tn C3N-00169_T  -l  druggability.log  -o alterations.out
    ```
    ```
    $ ./druggability.py  -t maf  -f LUAD.Somatic.050919.mnp.annot.maf -nn C3N-00169_N -tn C3N-00169_T  -l druggability.log  -o alterations.out  -at chol  -ato  trials.aux
    ```

2. **Ad hoc MAFs**. Variant call files, such as those redistributed with publications, may have a different format than those listed above. In that case, `druggability` can read a tab-delimited file whose first 13 fields are the following:
    ```
    Hugo_Symbol	NCBI_Build	Chromosome	Start_Position	End_Position	Variant_Classification	Variant_Type	Reference_Allele	Tumor_Seq_Allele1	Tumor_Sample_Barcode	Matched_Norm_Sample_Barcode	HGVSc	HGVSp_Short
    ```

    and analyzed with the `basicmaf` option:

    ```
    $ ./druggability.py  -t basicmaf  -f cohort.maf -nn 111 -tn 111  -l druggability.log  -o alterations.out  -at chol  -ato  trials.aux
    ```
    Notes:
	- Coordinates are 1-based
	- `Variant_Type` should be SNP, DNP, MNP, INS, or DEL

### Fusion example:
- Input file must be in ding-lab format.
    ```
    $ ./druggability.py  -t fusion -f CPTAC_fusions_v0.1.csv.tsv -tn C3L-00165_T  -l druggability.log  -o alterations.out
    ```
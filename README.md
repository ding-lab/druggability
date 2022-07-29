# druggability

## Requirements
- python3 environment
- alteration database
- variant file (single, or multisample, MAF or fusion results file; VCF not yet supported)
- (optional) clinical trials metadata file

## Installation
- `git clone https://github.com/ding-lab/druggability.git`
- `git clone https://github.com/ding-lab/druggability_databases.git`

## Configuration
In the [druggability_databases](https://github.com/ding-lab/druggability_databases/) repository, edit `config.py` with the correct paths. As the clinical trials metadata file is being updated and tested continually, it is being distributed through the repo at this time.

## Running
```
$ ./druggability.py -h
usage: druggability.py [-h] [-o OUTPUT_FILE] [-l LOG_FILE] [-d] [-nn NORMAL_NAME] [-at ANNOTATE_TRIALS]
                       [-ato TRIALS_AUXILIARY_OUTPUT_FILE] -t VARIATION_TYPE -f VARIANT_FILE -tn TUMOR_NAME

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
  -t VARIATION_TYPE     variation type: maf | fusion
  -f VARIANT_FILE       variant filename
  -tn TUMOR_NAME        sample or tumor sample name
```
Notes:
- The code is developed with cohort input files in mind and so it's best run as part of a loop. See the `examples` directory for script templates.
- For fusion inputs, the "tumor name" is thought of as any sample, normal or tumor.

### MAF examples:
- MAFs must be in ding-lab or PanCan union maf formats.
```
$ ./druggability.py  -t maf  -f LUAD.Somatic.050919.mnp.annot.maf -nn C3N-00169_N -tn C3N-00169_T  -l  druggability.log  -o alterations.out
```
```
$ ./druggability.py  -t maf  -f LUAD.Somatic.050919.mnp.annot.maf -nn C3N-00169_N -tn C3N-00169_T  -l druggability.log  -o alterations.out  -at chol  -ato  trials.aux
```

### Fusion example:
- Input file must be in ding-lab format.
```
$ ./druggability.py  -t fusion -f CPTAC_fusions_v0.1.csv.tsv -tn C3L-00165_T  -l druggability.log  -o alterations.out
```

# druggability

## Requirements
- python3 environment
- Installed repositories
- variant file (i.e., MAF)

## Installation
- `git clone https://github.com/ding-lab/druggability.git`
- `git clone https://github.com/ding-lab/druggability_databases.git`

## Configuration
In the `druggability` repo, edit `config.py` with the correct path to the druggability databases.

## Running
```
$ ./druggability.py -h
usage: druggability.py [-h] -t VARIATION_TYPE -f VARIANT_FILE

optional arguments:
  -h, --help         show this help message and exit

required arguments:
  -t VARIATION_TYPE  variation type: maf | fusion
  -f VARIANT_FILE    variant filename
```
Example:
```
./druggability.py  -t maf  -f LUAD.Somatic.050919.mnp.annot.maf   > out
```
**NOTE:** Fusions are a work-in-progress.

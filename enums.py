
# R. Jay Mashl <rmashl@wustl.edu>

# fusion filetypes
SINGLE     = 1    # single-sample report
COMBINED   = 2    # aggregated report

# maf filetypes
WASHU_MAF  = 1
UNION_MAF  = 2
BASIC_MAF  = 3

# variant classes
WILDTYPE   = 0
MUTATION   = 1
FUSION     = 2
INSERTION  = 3
DELETION   = 4
VARIANT_CLASSES = [WILDTYPE, MUTATION, FUSION, INSERTION, DELETION]

# trial alteration events
QUALIFYING    = 10
DISQUALIFYING = 11   # this value needs to be distinct from those in variant classes

# global
UNDECLARED = 0

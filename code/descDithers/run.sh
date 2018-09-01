#!/bin/bash

# source sims_maf
source /global/common/software/lsst/cori-haswell-gcc/stack/setup_current_sims.sh
# set up my own version rn
setup sims_maf -r /global/homes/a/awan/LSST/lsstRepos/sims_maf

# run the script
python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/code/descDithers/descDiths_wp_cadences.py
for file in `ls descDither*.csv`; do echo "gzipping $file"; gzip $file; done

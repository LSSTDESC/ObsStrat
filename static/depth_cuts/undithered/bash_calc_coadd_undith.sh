#!/bin/bash

# source sims_maf
source /global/common/software/lsst/cori-haswell-gcc/stack/setup_current_sims.sh
# set up my own version rn
setup sims_maf -r /global/homes/a/awan/LSST/lsstRepos/sims_maf
setup sims_maf_contrib -r /global/homes/a/awan/LSST/lsstRepos/sims_maf_contrib

##############################################################################
# run the script for all the cadences
for yr_cut in 1 3 6 10
do
    python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/code/modifyWFD/calc_coadd_depth.py \
                            --nside=256 --yr_cut=$yr_cut --no_dith --bands='i' \
                            --dbs_path='/global/cscratch1/sd/awan/dbs_wp_unzipped' \
                            --outDir='/global/cscratch1/sd/awan/lsst_output/coadd_output_noDith/' &
done
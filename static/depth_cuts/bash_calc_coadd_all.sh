#!/bin/bash

# source sims_maf
source /global/common/software/lsst/cori-haswell-gcc/stack/setup_current_sims.sh
# set up my own version rn
setup sims_maf -r /global/homes/a/awan/LSST/lsstRepos/sims_maf
setup sims_maf_contrib -r /global/homes/a/awan/LSST/lsstRepos/sims_maf_contrib

# run the script for all the cadences
#for yr_cut in 1 3 6 10
#do
#    python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/code/modifyWFD/calc_coadd_depth.py --nside=256 --yr_cut=$yr_cut \
#                                                            --outDir='/global/homes/a/awan/LSST/output/coadd_output_allwps_perNight/' &
#done

# run the script for the new cadences
for yr_cut in 1 3 6 10
do
    for db in kraken_2042.db kraken_2044.db mothra_2049.db nexus_2097.db
    do
        python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/code/modifyWFD/calc_coadd_depth.py --nside=256 --yr_cut=$yr_cut \
                                                            --specific_db=$db \
                                                            --outDir='/global/homes/a/awan/LSST/output/coadd_output_allwps_perNight/' &
    done
done
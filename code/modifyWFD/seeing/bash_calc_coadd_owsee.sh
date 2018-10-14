#!/bin/bash

# source sims_maf
source /global/common/software/lsst/cori-haswell-gcc/stack/setup_current_sims.sh
# set up my own version rn
setup sims_maf -r /global/homes/a/awan/LSST/lsstRepos/sims_maf
setup sims_maf_contrib -r /global/homes/a/awan/LSST/lsstRepos/sims_maf_contrib

# run the script for baseline and pontus2002 only
# dithered: 10yr
#for tag  in ow6 ow6c ow7 ow7c
#do
#python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/code/modifyWFD/calc_coadd_depth.py --nside=256 \
#                                                    --baseline_and_wide_only \
#                                                    --dbs_path=/global/cscratch1/sd/awan/owsee_dbs/${tag} \
#                                                    --outDir='/global/cscratch1/sd/awan/lsst_output/coadd_output_owsee/' &
#done

# run the script for Y1,3,6,10 for all cadences for ow6 only.
for tag  in ow6
do
    for yr_cut in 1 3 6 10
    do
        python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/code/modifyWFD/calc_coadd_depth.py --nside=256 \
                                                    --yr_cut=$yr_cut \
                                                    --dbs_path=/global/cscratch1/sd/awan/owsee_dbs/${tag} \
                                                    --outDir='/global/cscratch1/sd/awan/lsst_output/coadd_output_owsee/' &
    done
done
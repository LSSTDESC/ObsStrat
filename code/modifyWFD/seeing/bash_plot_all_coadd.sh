#!/bin/bash

# source sims_maf
source /global/common/software/lsst/cori-haswell-gcc/stack/setup_current_sims.sh
# set up my own version rn
setup sims_maf -r /global/homes/a/awan/LSST/lsstRepos/sims_maf

# plot the coadd maps for dbs with OpSims seeing
python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/code/modifyWFD/seeing/plot_coadd_maps.py \
                    --data_dir='/global/cscratch1/sd/awan/lsst_output/coadd_output_allwps_perNight' \
                    --outDir='/global/cscratch1/sd/awan/lsst_output/owsee_comparisons' &

# plot the coadd maps for Eric N.'s dbs
for tag  in ow6 ow6c ow7 ow7c
do
    python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/code/modifyWFD/seeing/plot_coadd_maps.py \
                        --data_dir='/global/cscratch1/sd/awan/lsst_output/coadd_output_owsee' \
                        --outDir='/global/cscratch1/sd/awan/lsst_output/owsee_comparisons' \
                        --owsee_tag=$tag &
done
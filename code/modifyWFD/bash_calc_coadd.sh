#!/bin/bash

# source sims_maf
source /global/common/software/lsst/cori-haswell-gcc/stack/setup_current_sims.sh
# set up my own version rn
setup sims_maf -r /global/homes/a/awan/LSST/lsstRepos/sims_maf
setup sims_maf_contrib -r /global/homes/a/awan/LSST/lsstRepos/sims_maf_contrib

# run the script
# dithered: 1yr
python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/code/modifyWFD/calc_coadd_depth.py --nside=256 \
                                                    --yr_cut=1 --baseline_and_wide_only \
                                                    --outDir='/global/homes/a/awan/LSST/output/coadd_output_allwps_perNight/' &
# dithered: 10yr
python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/code/modifyWFD/calc_coadd_depth.py --nside=256 \
                                                    --baseline_and_wide_only \
                                                    --outDir='/global/homes/a/awan/LSST/output/coadd_output_allwps_perNight/' &
# undithered : 1yr                                                             
python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/code/modifyWFD/calc_coadd_depth.py --nside=256 \
                                                    --yr_cut=1 --baseline_and_wide_only --no_dither \
                                                    --outDir='/global/homes/a/awan/LSST/output/coadd_output_noDith/' &
# undithered : 1yr                                                             
python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/code/modifyWFD/calc_coadd_depth.py --nside=256 \
                                                    --baseline_and_wide_only --no_dither \
                                                    --outDir='/global/homes/a/awan/LSST/output/coadd_output_noDith/' &
#!/bin/bash

########################################################################################################################
# The script implements depth + ebv cuts on ow6 coadd data for Y1,3,6,10.
########################################################################################################################
# source sims_maf
source /global/common/software/lsst/cori-haswell-gcc/stack/setup_current_sims.sh
# set up my own version rn
setup sims_maf -r /global/homes/a/awan/LSST/lsstRepos/sims_maf

##########################################################################################################
# run the script for all the cadences
for db in baseline2018a kraken_2026 kraken_2035 colossus_2665 \
          colossus_2664 colossus_2667 pontus_2002 pontus_2489 pontus_2502 mothra_2045
          #alt_sched alt_sched_rolling
do
    echo Working on $db
    python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/static/depth_cuts/implement_depth_ebv_cuts.py  \
                --db_path='/global/cscratch1/sd/awan/owsee_dbs/ow6/ow6_'${db}'.db' \
                --coadd_data_dir='/global/cscratch1/sd/awan/lsst_output/coadd_output_owsee/' \
                --yr_cuts='1yr, 3yr, 6yr, 10yr' \
                --chosen_cuts='24.5, 25.0, 25.5, 26.0' \
                --ebv_cut --save_stuff --dont_show_plots \
                --outDir='/global/cscratch1/sd/awan/lsst_output/owsee_comparisons/depth_cuts/' \
                --outDir_md='/global/homes/a/awan/LSST/lsstRepos/ObsStrat/code/modifyWFD/seeing/mdfiles/'
done
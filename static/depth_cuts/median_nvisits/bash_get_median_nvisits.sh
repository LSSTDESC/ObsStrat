#!/bin/bash

# source sims_maf
source /global/common/software/lsst/cori-haswell-gcc/stack/setup_current_sims.sh
# set up my own version rn
setup sims_maf -r /global/homes/a/awan/LSST/lsstRepos/sims_maf

##############################################################################
# run the script for all the cadences
for db in baseline2018a kraken_2026 kraken_2035 kraken_2036 colossus_2665 \
          colossus_2664 colossus_2667 pontus_2002 pontus_2489 pontus_2502 mothra_2045 \
           kraken_2042 kraken_2044 mothra_2049 nexus_2097
do
    echo Working on $db
    python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/static/depth_cuts/median_nvisits/get_median_nvisits.py  \
                --db_path='/global/cscratch1/sd/awan/dbs_wp_unzipped/'${db}'.db' \
                --outDir='/global/homes/a/awan/desc/depth_data_outputs/' \
                --eg_path='/global/homes/a/awan/desc/depth_data_outputs/' &
done

#for db in alt_sched alt_sched_rolling
#do
#    echo Working on $db
#    python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/static/depth_cuts/median_nvisits/get_median_nvisits.py  \
#                --altsched \
#                --db_path='/global/cscratch1/sd/awan/dbs_wp_unzipped/slair_altsched/'${db}'.db' \
#                --outDir='/global/homes/a/awan/desc/depth_data_outputs/' \
#                --eg_path='/global/homes/a/awan/desc/depth_data_outputs/' &
#done
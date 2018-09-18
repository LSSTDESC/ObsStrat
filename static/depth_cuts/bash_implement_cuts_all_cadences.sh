#!/bin/bash

# source sims_maf
source /global/common/software/lsst/cori-haswell-gcc/stack/setup_current_sims.sh
# set up my own version rn
setup sims_maf -r /global/homes/a/awan/LSST/lsstRepos/sims_maf

current_dir=$(pwd)
# run the script for all the cadences
#for db in baseline2018a kraken_2026 kraken_2035 kraken_2036 colossus_2665 \
#          colossus_2664 colossus_2667 pontus_2002 pontus_2489 pontus_2502 mothra_2045
for db in kraken_2042 kraken_2044 mothra_2049 nexus_2097
do
    echo Working on $db
    python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/static/depth_cuts/implement_depth_ebv_cuts.py  \
                --db_path='/global/cscratch1/sd/awan/dbs_wp_unzipped/'${db}'.db' \
                --coadd_data_dir='/global/homes/a/awan/LSST/output/coadd_output_allwps_perNight/' \
                --yr_cuts='1yr, 3yr, 6yr, 10yr' \
                --chosen_cuts='24.5, 25.0, 25.5, 26.0' \
                --ebv_cut --save_stuff --dont_show_plots \
                --outDir='/global/homes/a/awan/desc/depth_data_outputs/' \
                --outDir_md='/global/homes/a/awan/LSST/lsstRepos/ObsStrat/static/depth_cuts/mdfiles/'
done

# change permissions on the outputs 
cd /global/homes/a/awan/desc/
chgrp -R lsst depth_data_outputs
chmod -R g-w depth_data_outputs
cd $current_dir
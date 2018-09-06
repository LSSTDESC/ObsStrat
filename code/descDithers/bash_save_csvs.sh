#!/bin/bash

# source sims_maf
source /global/common/software/lsst/cori-haswell-gcc/stack/setup_current_sims.sh
# set up my own version rn
setup sims_maf -r /global/homes/a/awan/LSST/lsstRepos/sims_maf

# run the script
python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/code/descDithers/save_csv_dithers.py \
                            --dbs_path='/global/cscratch1/sd/awan/dbs_wp_unzipped' \
                            --outDir='/global/homes/a/awan/desc/wp_descDithers_csvs' \
                            --rot_rand_seed=42 --trans_rand_seed=42 --save_plots
                            #--db_files_only='baseline2018a.db'
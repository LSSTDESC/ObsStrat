#!/bin/bash
##############################################################################################################
# The goal here is to run save_plot_rot_data_newfbs for various dbs to save rotational angle related data.
#
# Humna Awan: humna.awan@rutgers.edu
#
##############################################################################################################
# source sims_maf
source /global/common/software/lsst/cori-haswell-gcc/stack/setup_current_sims.sh
setup lsst_sims
# set up my mafcontrib
setup sims_maf_contrib -r /global/homes/a/awan/LSST/lsstRepos/sims_maf_contrib

########################################################################################################
dbs_path=/global/cscratch1/sd/awan/dbs_post_wp/
#######################################################################################
# run the script for the various OpSim outputs for non-DD visits.
outdir=/global/cscratch1/sd/awan/lsst_output/post_wp_output/rot_output/
for db_path in $(find ${dbs_path} -name '*.db' )
do
    echo 'Saving rot data for '${db_path}
    python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/code/descDithers/save_plot_rot_data_newfbs.py \
                        --outdir=${outdir} --dbfile=${db_path}
done

#######################################################################################
outdir=/global/cscratch1/sd/awan/lsst_output/post_wp_output/rot_output_dd/
# run for just the DD visits
for db_path in $(find ${dbs_path} -name '*.db' )
do
    echo 'Saving rot data for '${db_path}
    python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/code/descDithers/save_plot_rot_data_newfbs.py \
                        --outdir=${outdir} --dbfile=${db_path} --dd
done
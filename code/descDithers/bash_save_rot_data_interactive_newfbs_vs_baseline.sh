#!/bin/bash
#####################################################################################################
# The goal here is to run save_rot_data.py and save_plot_rot_data_newfbs.py to save rot-data for
# the baseline and and the FBS output with rotational dithers.
# Humna Awan: humna.awan@rutgers.edu
#
#####################################################################################################
# source sims_maf from the desc stack
source /global/common/software/lsst/cori-haswell-gcc/stack/setup_current_sims.sh
#setup lsst_sims
# set up my maf
setup sims_maf -r /global/homes/a/awan/LSST/lsstRepos/sims_maf
#setup sims_maf_contrib -r /global/homes/a/awan/LSST/lsstRepos/sims_maf_contrib

outdir=/global/homes/a/awan/desc/rot_output/

#######################################################################################
# run the script for the various OpSim outputs.
for db in baseline2018a kraken_2026
do
    # script for saving the data
    python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/code/descDithers/save_rot_data.py \
                        --outdir=${outdir} \
                        --dbfile='/global/cscratch1/sd/awan/dbs_wp_unzipped/'${db}'.db'

    # script for plotting the saved data
    python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/code/descDithers/save_rot_plot.py \
                        --outdir=${outdir} \
                        --dbname=${db}
done

#######################################################################################
# run the script for the new FBS output
python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/code/descDithers/save_plot_rot_data_newfbs.py \
                        --outdir=${outdir} \
                        --dbfile='/global/cscratch1/sd/awan/dbs_wp_unzipped/rotator_1exp_pairsmix_10yrs.db'
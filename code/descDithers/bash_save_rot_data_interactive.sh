#!/bin/bash
#####################################################################################################
# The goal here is to run save_rot_data.py to save data for the dithered/undithed rotational angles
# for different dbs. 
#
# Not considering alt_sched and FBS outputs since they don't have fieldIDs (and FBS dont seem to have
# the LST column needed the PA stacker).
#
# Humna Awan: humna.awan@rutgers.edu
#
#####################################################################################################
# source sims_maf from the desc stack
source /global/common/software/lsst/cori-haswell-gcc/stack/setup_current_sims.sh
#setup lsst_sims
# set up my maf
setup sims_maf -r /global/homes/a/awan/LSST/lsstRepos/sims_maf
#setup sims_maf_contrib -r /global/homes/a/awan/LSST/lsstRepos/sims_maf_contrib

outdir=/global/homes/a/awan/desc/wp_descDithers_csvs/compare_rot_dith/

#######################################################################################
# run the script for the various OpSim outputs.
for db in baseline2018a pontus_2002 \
            kraken_2026 kraken_2035 colossus_2665 \
            colossus_2664 colossus_2667 pontus_2489 \
            kraken_2042 kraken_2044 mothra_2049 nexus_2097\
            pontus_2502 kraken_2036 mothra_2045
do
    # script for saving the data
    #python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/code/descDithers/save_rot_data.py \
    #                    --outdir=${outdir} \
    #                    --dbfile='/global/cscratch1/sd/awan/dbs_wp_unzipped/'${db}'.db' &

    # script for plotting the saved data
    python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/code/descDithers/save_rot_plot.py \
                        --outdir=${outdir} \
                        --dbname=${db} &
done
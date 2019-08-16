#!/bin/bash
##############################################################################################################
# The goal here is to run the fom_emulator for various dbs.
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
fom_path=/global/homes/a/awan/LSST/lsstRepos/ObsStrat/static/fom_emulator/FoM
data_file_tag=stats
renorm_strategy=baseline10yrs
# all dbs
data_path=/global/homes/a/awan/LSST/lsstRepos/ObsStrat/postwp/results/stats/
fig_dir=/global/homes/a/awan/LSST/lsstRepos/ObsStrat/postwp/results/fom_figs/

# no constraints on the color bar, etc for the contour plots.
python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/static/fom_emulator/fom_emulator.py \
                --fom_path=${fom_path} --data_path=${data_path} \
                --data_file_tag=${data_file_tag} --fig_dir=${fig_dir} \
                --renorm_strategy=${renorm_strategy} --group_dbs

# impplement constraints on the color bar, etc for the contour plots.
python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/static/fom_emulator/fom_emulator.py \
                --fom_path=${fom_path} --data_path=${data_path} \
                --data_file_tag=${data_file_tag} --fig_dir=${fig_dir} \
                --renorm_strategy=${renorm_strategy} --group_dbs --limit_contour
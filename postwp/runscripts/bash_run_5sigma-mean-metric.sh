#!/bin/bash
##############################################################################################################
# The goal here is to run calc_mena_col for 5sigma for various dbs.
#
# Humna Awan: humna.awan@rutgers.edu
#
##############################################################################################################
# source sims_maf
source /global/common/software/lsst/cori-haswell-gcc/stack/setup_current_sims.sh
setup lsst_sims
# set up my mafcontrib
# setup sims_maf -r /global/homes/a/awan/LSST/lsstRepos/sims_maf
setup sims_maf_contrib -r /global/homes/a/awan/LSST/lsstRepos/sims_maf_contrib

########################################################################################################
outdir=/global/cscratch1/sd/awan/lsst_output/post-wp_v1.4_five-sigma-depth/
dbs_path=/global/cscratch1/sd/awan/dbs_post_wp_v3/

for db_path in $(find ${dbs_path} -name '*.db')
do
    echo 'Running analysis for '${db_path}
    python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/postwp/runscripts/calc_mean_col.py \
                                --nside=256 \
                                --db_path=${db_path} \
                                --outdir=${outdir} \
                                --colname='fiveSigmaDepth' \
                                --band='i'
    done
done

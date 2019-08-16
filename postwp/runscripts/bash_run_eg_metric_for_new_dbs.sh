#!/bin/bash
##############################################################################################################
# The goal here is to run run_eg_metric for various dbs.
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
outdir=/global/cscratch1/sd/awan/lsst_output/post_wp_output/
dbs_path=/global/cscratch1/sd/awan/dbs_post_wp/

for yr_cut in 10 1 3 6
do
    for db_path in $(find ${dbs_path} -name '*.db')
    do
        echo 'Running analysis for '${yr_cut}'yr for '${db_path}
        python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/postwp/runscripts/run_eg_metric.py \
                            --outdir=${outdir} \
                            --db_path=${db_path} \
                            --nside=256 --yr_cut=${yr_cut} --fbs
    done
done


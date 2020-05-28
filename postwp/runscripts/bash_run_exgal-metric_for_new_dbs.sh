#!/bin/bash
##############################################################################################################
# The goal here is to run run_eg_metric for various dbs.
#
# Humna Awan: humna.awan@rutgers.edu
#
##############################################################################################################
# source sims_maf
#source /global/common/software/lsst/cori-haswell-gcc/stack/setup_current_sims.sh
#setup lsst_sims
# set up my mafcontrib
# setup sims_maf -r /global/homes/a/awan/LSST/lsstRepos/sims_maf
#setup sims_maf_contrib -r /global/homes/a/awan/LSST/lsstRepos/sims_maf_contrib

########################################################################################################
maindir=/global/cscratch1/sd/awan/lsst_output/post_wp_output_v1.4_-0.1cuts
outdir=${maindir}/exgalm5_data/
#outdir=/global/homes/a/awan/LSST/lsstRepos/ObsStrat/postwp/runscripts/test_exgalm5_data/
dbs_path=/global/cscratch1/sd/awan/dbs_post_wp_v3/
egpath=${maindir}/lss_fom/eg-footprint-mask/

for yr_cut in 10 #1 3 6
do
    for band in i u g r z y
    do
        for db_path in $(find ${dbs_path} -name '*.db')
        do
            echo 'Running analysis for '${band}'-band for '${yr_cut}'yr for '${db_path}
            python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/postwp/runscripts/run_exgal_metric.py \
                                --outdir=${outdir} \
                                --eg-path=${egpath} \
                                --db-path=${db_path} \
                                --nside=256 --yr_cut=${yr_cut} --band=${band}
        done
    done
done

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
setup sims_maf_contrib -r /global/homes/a/awan/LSST/lsstRepos/sims_maf_contrib

########################################################################################################
outdir=/global/cscratch1/sd/awan/lsst_output/post_wp_output/older_dbs/

for db in baseline2018a pontus_2002 kraken_2026 colossus_2665 mothra_2045
do
    for yr_cut in 1 10
    do
        echo 'Running analysis for '${yr_cut}'yr for '${db}
        python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/code/postwp/runscripts/run_eg_metric.py \
                            --outdir=${outdir} \
                            --db_path='/global/cscratch1/sd/awan/dbs_wp_unzipped/'${db}'.db' \
                            --nside=128 --yr_cut=${yr_cut} --wdith --older_cuts
    done
done

# now run the y1 data with new depth cut
for db in baseline2018a pontus_2002 kraken_2026 colossus_2665 mothra_2045
do
    for yr_cut in 1
    do
        echo 'Running analysis for '${yr_cut}'yr for '${db}
        python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/code/postwp/runscripts/run_eg_metric.py \
                            --outdir=${outdir} \
                            --db_path='/global/cscratch1/sd/awan/dbs_wp_unzipped/'${db}'.db' \
                            --nside=128 --yr_cut=${yr_cut} --wdith
    done
done


#!/bin/bash
##############################################################################################################
# The goal here is to run get_overlaps for various dbs.
########################################################################################################
fbs_version=v1.7
maindir=/global/cscratch1/sd/awan/lsst_output/post_wp_output_${fbs_version}_-0.1cuts-1-olderAx/
repo_path=/global/homes/a/awan/LSST/lsstRepos/ObsStrat/postwp
# path for the dbs order list
dbs_order_path=${repo_path}/paper-data/summary_csv_${fbs_version}_-0.1cuts/given_order_${fbs_version}.csv
# lsst footprint details
eg_path=${maindir}/eg-footprint-mask/
nside_lsst=256
# outdir
outdir=${repo_path}/paper-data/overlaps_${fbs_version}
# ------------------------------------------------------------------------
# run things for 4MOST
survey_name=4MOST
other_footprint_path=/global/homes/a/awan/desc_oswg/2021-4most_nside16.csv
nside_other=16

echo 'Running analysis for '${fbs_version} 'dbs for overlaps with '${survey_name}
python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/postwp/runscripts/get_overlaps.py \
                        --outdir=${outdir} \
                        --eg-path=${eg_path} --nside-lsst=${nside_lsst} \
                        --other-footprint-path=${other_footprint_path} \
                        --survey=${survey_name} --nside-other=${nside_other} \
                        --dbs-order-path=${dbs_order_path}
# ------------------------------------------------------------------------
# run things for DESI
survey_name=DESI
other_footprint_path=/global/homes/a/awan/desc_oswg/4most-tides_desi_data/DESI_pixels_nside256_ring.csv
nside_other=256

echo 'Running analysis for '${fbs_version} 'dbs for overlaps with '${survey_name}
python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/postwp/runscripts/get_overlaps.py \
                        --outdir=${outdir} \
                        --eg-path=${eg_path} --nside-lsst=${nside_lsst} \
                        --other-footprint-path=${other_footprint_path} \
                        --survey=${survey_name} --nside-other=${nside_other} \
                        --dbs-order-path=${dbs_order_path}
# ------------------------------------------------------------------------
# run things for euclid
survey_name=euclid
other_footprint_path=/global/homes/a/awan/desc_oswg/euclid_footprint_pixels_nside256.csv
nside_other=256

echo 'Running analysis for '${fbs_version} 'dbs for overlaps with '${survey_name}
python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/postwp/runscripts/get_overlaps.py \
                        --outdir=${outdir} \
                        --eg-path=${eg_path} --nside-lsst=${nside_lsst} \
                        --other-footprint-path=${other_footprint_path} \
                        --survey=${survey_name} --nside-other=${nside_other} \
                        --dbs-order-path=${dbs_order_path}
# ------------------------------------------------------------------------
# create plot
python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/postwp/runscripts/plot_overlaps.py \
                        --fbs-version=${fbs_version}

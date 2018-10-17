#!/bin/bash
########################################################################################################################
# The goal here is to check whether we have the files for different cadences, for different
# years for different bands. If the coadd file isn't found, a message is printed.
#
# This script is helpful since the coadd runs (by bash_calc_coadd) don't always complete successfully.
# It checks the availability of needed files before running bash_implement_cuts which needs coadd outputs.
########################################################################################################################
coadd_data_dir=/global/cscratch1/sd/awan/lsst_output/coadd_output_allwps_slair_altsched

for yr_cut in 1yearCut 3yearCut 6yearCut fullSurveyPeriod
do
    # loop over bands
    for band in u g r i z y
    do
        # loop over different cadences
        for db in cadence_roll_75_mix_rolling_mix_10yrs \
                    roll_mix_100_rolling_mix_10yrs \
                    roll_mix_rolling_mix_10yrs \
                    rolling_10yrs  tms_roll_10yrs \
                    alt_sched  alt_sched_rolling
        do
            path=coaddM5Analysis_nside256_withDustExtinction_0pixelRadiusForMasking_${band}Band_${db}_${yr_cut}_directory
            #echo ${coadd_data_dir}/${path}/unmaskedCoaddData/*.npz
            count=`ls -1 ${coadd_data_dir}/${path}/unmaskedCoaddData/*.npz 2>/dev/null | wc -l`
            #echo $count
            if [ $count == 0  ]; then
                echo $count
                echo "No npz file not found!"
                echo ${coadd_data_dir}/${path}/unmaskedCoaddData/*.npz
            fi
        done
    done
done
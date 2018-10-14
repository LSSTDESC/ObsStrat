#!/bin/bash

cd /global/cscratch1/sd/awan/dbs_wp_unzipped/slair_altsched/

##############################################################################
# get slair outputs
for i in rolling/rolling_10yrs tms_roll/tms_roll_10yrs
do
	echo getting $i
	wget https://lsst-web.ncsa.illinois.edu/sim-data/beta_slair_surveys/runs/$i.db
done


# get the rolling mix ones. need to rename the db.
for i in cadence_roll_75_mix roll_mix_100 roll_mix
do
	echo getting $i/rolling_mix_10yrs.db
	wget https://lsst-web.ncsa.illinois.edu/sim-data/beta_slair_surveys/runs/$i/rolling_mix_10yrs.db
    mv rolling_mix_10yrs.db "$i"_rolling_mix_10yrs.db
done

##############################################################################
# get alt_sched outputs
for i in alt_sched alt_sched_rolling
do
	echo getting $i
	wget http://altsched.rothchild.me:8080/$i.db
done

##############################################################################
# change permissions
cd /global/cscratch1/sd/awan/dbs_wp_unzipped/
chgrp -R lsst slair_altsched/
chmod -R g-w slair_altsched/

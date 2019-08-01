#!/bin/bash

# colossus_2665 pontus_2002 colossus_2664 colossus_2667 pontus_2489 kraken_2035 mothra_2045 pontus_2502 kraken_2036 kraken2026 baseline2018a
for i in kraken_2042 kraken_2044 mothra_2049 nexus_2097
do
	echo starting with: $i
	cd /global/cscratch1/sd/awan/dbs_wp_unzipped/
	wget http://astro-lsst-01.astro.washington.edu:8080/db_gzip/$i.db.gz

	gunzip $i.db.gz
	chgrp lsst $i.db
	chmod g-w $i.db
done




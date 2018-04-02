#!/bin/sh
# This shell script converts the text (tab separated values) table
# produted by simsee to an sqlite database that can be used by
# opsim4

INFILE=${1}
OUTFILE=${2}
CSVFILE=$(mktemp)
tail -n +2 ${INFILE} | cat -n | cut -f1,3,5  --output-delimiter="," > ${CSVFILE}

sqlite3 $OUTFILE <<EOF
CREATE TABLE Seeing(seeingId INTEGER PRIMARY KEY,s_date INTEGER,seeing DOUBLE);
.mode csv
.import ${CSVFILE} Seeing
EOF

rm $CSVFILE

##############################################################################################################
# This scripts saves dithered and undithered rotTelPos and rotSkyPos for different dbs. The code was initially
# run in an iPython notebook but it takes a really long time to load the dbs; its better to just save the data
# in one go to allow quick reproducibility.
#
# Humna Awan: humna.awan@rutgers.edu
#
##############################################################################################################
import os
import pandas as pd
import lsst.sims.maf.db as db
import lsst.sims.maf.stackers as stackers
import time
#######################################################################################
from optparse import OptionParser
parser = OptionParser()
parser.add_option('--dbfile', dest='dbfile',
                  help='Path to the folder with the db to consider.',
                  default='/global/cscratch1/sd/awan/dbs_wp_unzipped/baseline2018a.db')
parser.add_option('--outdir', dest='outdir',
                  help='Path to the folder with the db to consider.',
                  default='/global/homes/a/awan/desc/wp_descDithers_csvs/compare_rot_dith/')

(options, args) = parser.parse_args()
dbfile = options.dbfile
outdir = options.outdir
#######################################################################################
rot_rand_seed = 42

dbname = dbfile.split('/')[-1].split('.db')[0]
print('Running save_rot_data for %s'%(dbname))
time0 = time.time()
# ----------------------------------------------------------
# connect to the database
opsdb = db.OpsimDatabase(dbfile)
# ----------------------------------------------------------
# WFD only
prop_ids, prop_tags = opsdb.fetchPropInfo()
wfd_constraint = opsdb.createSQLWhere('WFD', prop_tags)

# fetch the data: columns need for the rotational dither and parallactic angle stacker
colnames=['fieldRA', 'fieldDec', 'fieldId', 'observationId', \
          'rotTelPos', 'rotSkyPos', 'night', 'filter', \
          'observationStartMJD', 'observationStartLST']
simdata = opsdb.fetchMetricData(colnames=colnames, sqlconstraint=wfd_constraint)

# ----------------------------------------------------------
# add rotational dithers: adds randomDitherPerFilterChangeRotTelPos column to simdata
s = stackers.RandomRotDitherPerFilterChangeStacker(degrees=opsdb.raDecInDeg,
                                                   randomSeed=rot_rand_seed)
simdata = s.run(simdata)

# ----------------------------------------------------------
# add parallactic angle values to the stacker: added PA column
# needed to go from dithered rotTelPos to dithered rotSkyPos; rotSkyPos = rotTelPos - PA
s = stackers.ParallacticAngleStacker(degrees=opsdb.raDecInDeg)
simdata = s.run(simdata)

# ----------------------------------------------------------
data_dict = {}
for col in ['observationId', 'fieldId', 'rotTelPos', \
            'randomDitherPerFilterChangeRotTelPos', \
            'rotSkyPos', 'PA']:
    data_dict[col] = simdata[col]

# ----------------------------------------------------------
# save the output array
filename = '%s_data.csv'%(dbname)
df = pd.DataFrame(data_dict)
df.to_csv('%s/%s'%(outdir, filename), index=False)
print('\nSaved %s'%filename)

print('Total time taken: %.2f min'%((time.time()-time0)/60.))
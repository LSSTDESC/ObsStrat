##############################################################################################################
# This scripts saves dithered and undithered rotTelPos and rotSkyPos for different dbs. The code was initially
# run in an iPython notebook but its takes really long to load the dbs; its better to just save the data
# for saving time and allowing quick reproducibility.
#
# Humna Awan: humna.awan@rutgers.edu
#
##############################################################################################################
import os
import pandas as pd
import lsst.sims.maf.db as db
import lsst.sims.maf.stackers as stackers
import time

dbs_path = '/global/cscratch1/sd/awan/dbs_wp_unzipped'
rot_rand_seed = 42
outDir = '/global/homes/a/awan/desc/wp_descDithers_csvs/compare_rot_dith/'

# get all the files
dbfiles = [f for f in os.listdir(dbs_path) if f.endswith('db')]

time0 = time.time()
# now read in the data
simdatas = {}
for i, dbfile in enumerate(dbfiles): # loop over all the db files
    print(dbfile)
    time1 = time.time()
    opsdb = db.OpsimDatabase('%s/%s'%(dbs_path, dbfile)) # connect to the database

    # WFD only
    prop_ids, prop_tags = opsdb.fetchPropInfo()
    wfd_constraint = opsdb.createSQLWhere('WFD', prop_tags)
    # fetch the data: columns need for the rotational dither and parallactic angle stacker
    colnames=['fieldRA', 'fieldDec', 'rotTelPos', 'rotSkyPos', 'night', 'filter',
              'observationStartMJD', 'observationStartLST']
    simdata = opsdb.fetchMetricData(colnames=colnames, sqlconstraint=wfd_constraint)
    
    # add rotational dithers: adds randomDitherPerFilterChangeRotTelPos column to simdata
    s = stackers.RandomRotDitherPerFilterChangeStacker(degrees=opsdb.raDecInDeg,
                                                       randomSeed=rot_rand_seed)
    simdata = s.run(simdata)
    
    # add parallactic angle values to the stacker: added PA column
    # needed to go from dithered rotTelPos to dithered rotSkyPos since rotSkyPos = rotTelPos - PA
    s = stackers.ParallacticAngleStacker(degrees=opsdb.raDecInDeg)
    simdata = s.run(simdata)
    # store the data
    key = dbfile.split('.db')[0]
    simdatas[key] = {}
    for col in ['rotTelPos', 'randomDitherPerFilterChangeRotTelPos', 'rotSkyPos', 'PA']:
        simdatas[key][col] = simdata[col]
    # save the output array
    filename = '%s_data.csv'%(key)
    df = pd.DataFrame(simdatas[key])
    df.to_csv('%s/%s'%(outDir, filename), index=False)
    print('\nSaved %s'%filename)
    
    print('Time taken for %s: %.2f min'%(dbfile, (time.time()-time1)/60.))

print('Total time taken: %.2f min'%((time.time()-time0)/60.))

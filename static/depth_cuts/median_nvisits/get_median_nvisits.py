# Goal here is to produce calculate the median number of visits
#
# Humna Awan: humna.awan@rutgers.edu
#
########################################################################
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

import lsst.sims.maf.db as db
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.metricBundles as metricBundles
import lsst.sims.maf.stackers as mafStackers
import time
import healpy as hp
import numpy as np
import pandas as pd

########################################################################
from optparse import OptionParser
parser = OptionParser()
parser.add_option('--db_path', dest='db_path',
                  help='Path to OpSim db to consider.',
                  default='/global/cscratch1/sd/awan/dbs_wp_unzipped/baseline2018a.db')
parser.add_option('--outDir', dest='outDir',
                  help='Path to the folder where all the output should be stored. The folder must already exist.',
                  default='/global/cscratch1/sd/awan/lsst_output/nvisits/')
parser.add_option('--eg_path', dest='eg_path',
                  help='Path to the folder where the depth/ebv-limited footprint is.',
                  default=None)
parser.add_option('--altsched',
                  action='store_true', dest='altsched', default=False,
                  help= 'Use the tag to specify a altscched run.')

(options, args) = parser.parse_args()
db_path = options.db_path
outDir = options.outDir
eg_path = options.eg_path
altsched = options.altsched
########################################################################
start_time = time.time()
dbname = db_path.split('/')[-1].split('.db')[0]
print('## Running the code for %s'%dbname)

# full footprint; no depth/ebv cuts
opsdb = db.OpsimDatabase(db_path)
# check if fieldId is present
if 'fieldId' not in opsdb.columnNames:
    raise ValueError('Dont have fieldIds for %s.'%(dbname))
# WFD only
propIds, propTags = opsdb.fetchPropInfo()
wfdWhere = opsdb.createSQLWhere('WFD', propTags)
# fetch the data
simdata = opsdb.fetchMetricData(colnames=['fieldId','fieldRA', 'night', 'fieldDec'], sqlconstraint=wfdWhere)
data = pd.DataFrame(simdata)

if len(np.unique(data['fieldId']))<=1:
    # something isnt right.
    raise ValueError('Somethings wrong: only found %s fieldIds for %s.'%(len(np.unique(data['fieldId'])), dbname))

filename = 'median_nvisits_per_field_10yr.csv'
if os.path.isfile('%s/%s'%(outDir, filename)):
    ofile = open('%s/%s'%(outDir, filename),'a')
else:
    ofile = open('%s/%s'%(outDir, filename),'w')
    ofile.write('dbname,median number of visits\n')
ofile.write('%s,%s\n'%(dbname, np.median(data.groupby('fieldId').size())))
ofile.close()

print('Wrote to %s'%filename)
print('## Time taken: %.2f min'%((time.time()-start_time)/60.))

# ----------------------------------------------------------------------------------------------------
if eg_path is not None:
    print('\n## Looking at the fieldIDs in the extragalactic footprint ... ')
    # read in the save pixels
    files = [f for f in os.listdir(eg_path) if f.endswith('csv') \
                                            and f.__contains__('_10yr_') \
                                            and f.__contains__('_%s_'%dbname)]
    print('%s files found: %s'%(len(files), files))
    file = files[0]
    data_eg = pd.read_csv('%s/%s'%(eg_path, file))
    nside = int(file.split('nside')[-1].split('_')[0])
    # now find all the fieldIDs that have observations
    # first need to set up the slicer
    if altsched:
        slicer = slicers.HealpixSlicer(lonCol='fieldRA',
                                       latCol='fieldDec',
                                       latLonDeg=opsdb.raDecInDeg, nside=nside)
    else:
        slicer = slicers.HealpixSlicer(lonCol='randomDitherPerNightRa',
                                       latCol='randomDitherPerNightDec',
                                       latLonDeg=opsdb.raDecInDeg, nside=nside)
        s = mafStackers.RandomDitherPerNightStacker(degrees=opsdb.raDecInDeg, randomSeed=1000)
        simdata = s.run(simdata)
    # slice the data now
    slicer.setupSlicer(simdata)
    # now find the fieldIds.
    fID_list = []
    for pixel in data_eg['pixNum']:
        indObsInPixel = slicer._sliceSimData(pixel)
        fID_list += list(simdata[indObsInPixel['idxs']]['fieldId']) # fieldIDs corresponding to pixel
    fID_list = np.unique(fID_list)

    all_fids = np.unique(data['fieldId'])
    print('Found %s fIDs in EG region vs. %s in total'%(len(fID_list), len(all_fids)))

    # now find the median across numbers for only the fieldIds in our list
    nvisits_all = np.array(data.groupby('fieldId').size())
    if len(all_fids)==len(fID_list):
        nvisit_subset = nvisits_all
    else:
        nvisit_subset = np.zeros(len(fID_list))
        for i, fID in enumerate(fID_list):
            ind = np.where(all_fids == fID)[0]
            nvisit_subset[i] = nvisits_all[ind]

    filename = 'median_nvisits_per_field_10yr_egfootprint.csv'
    if os.path.isfile('%s/%s'%(outDir, filename)):
        ofile = open('%s/%s'%(outDir, filename),'a')
    else:
        ofile = open('%s/%s'%(outDir, filename),'w')
        ofile.write('dbname,median number of visits\n')
    ofile.write('%s,%s\n'%(dbname, np.median(nvisit_subset)))
    ofile.close()

    print('Wrote to %s'%filename)
    print('## Total time taken: %.2f min\n'%((time.time()-start_time)/60.))
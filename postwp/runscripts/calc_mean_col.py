# the goal here is to calculate mean of the specified column in the specified db across non-dd Y10 visits.
# save the bundle data as an npz file in the specified outdir.
#
# Humna Awan: humna.awan@rutgers.edu
#
########################################################################################################################
import matplotlib
matplotlib.use('Agg')
import os
import time
import numpy as np

import lsst.sims.maf.metricBundles as metricBundles
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.db as db

########################################################################################################################
from optparse import OptionParser
parser = OptionParser()
parser.add_option('--nside', dest='nside', type='int',
			help='HEALPix resolution parameter.', default=256)
parser.add_option('--db_path', dest='db_path',
                  help='Path to the db file.')
parser.add_option('--outdir', dest='outdir',
                  help='Path to the folder where all the output should be stored. The folder must already exist.')
parser.add_option('--colname', dest='colname',
                  help='Column name for which to save data.')
parser.add_option('--band', dest='band',
                  help='Bands to consider.')

(options, args) = parser.parse_args()
nside = options.nside
db_path = options.db_path
outdir = options.outdir
band = options.band
colname = options.colname

yr_cut = 10
########################################################################################################################
print('## running calc_mean_cols.py for %s band for %s' % (band, colname) )
time0 = time.time()
# make the outdir
os.makedirs(outdir, exist_ok=True)

# get the dbname
dbname = db_path.split('/')[-1].split('.db')[0]
# connect to the database
opsdb = db.OpsimDatabase(db_path)
# ----------------------------------------------------------
# non-DD visits
# assemble the sql constraint
sqlconstraint = 'night <= %s'%(yr_cut * 365.25)
sqlconstraint += ' and note not like "DD%"'

# set up the slicer
slicer = slicers.HealpixSlicer(lonCol='fieldRA', latCol='fieldDec',
                               latLonDeg=opsdb.raDecInDeg, nside=nside)
# the resultsdb object
resultsDb = db.ResultsDb(database='resultsDb.db', outDir=outdir)
# set up the mean metric
mean_metric = metrics.MeanMetric(col=colname)
# set up the metric bundles object
bundle = metricBundles.MetricBundle(metric=mean_metric,
                                    slicer=slicer,
                                    sqlconstraint=sqlconstraint,
                                    runName=dbname)
grp = metricBundles.MetricBundleGroup(bundleDict={0: bundle}, dbObj=opsdb,
                                      outDir=outdir,
                                      resultsDb=resultsDb, saveEarly=False)
grp.runAll()

# ----------------------------------------------------------
# save the data array
filename = '%s_mean-%s_y%s_%sband_nside%s.npz'%(dbname, colname, yr_cut, band, nside)

bundle.slicer.writeData(outfilename='%s/%s'%(outdir, filename),
                        metricValues=bundle.metricValues,
                        metricName=bundle.metric.name,
                        simDataName=bundle.runName,
                        constraint=bundle.constraint,
                        metadata=bundle.metadata,
                        displayDict=bundle.displayDict,
                        plotDict=bundle.plotDict)
print('\n## saved %s' % filename)

print('## total time taken: %.2f min'%((time.time()-time0)/60.))
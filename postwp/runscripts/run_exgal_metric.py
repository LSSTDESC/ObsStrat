##############################################################################################################
# This script runs the the exgal metric for the specified db; the data is saved as a metricBundle.
# Then some statistics are calculate in the egfootprint, whose main path is user-specified.
#
# Humna Awan: humna.awan@rutgers.edu
#
##############################################################################################################
import numpy as np
import os
import healpy as hp
import time

import lsst.sims.maf.maps as maps
import lsst.sims.maf.db as db
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.metricBundles as metricBundles

#######################################################################################
from optparse import OptionParser
parser = OptionParser()
parser.add_option('--outdir', dest='outdir',
                  help='Path to the folder where all the output should be stored; should already exist.')
parser.add_option('--db_path', dest='db_path',
                  help='Path to the folder with the db to consider.')
parser.add_option('--eg-path', dest='eg_path',
                  help='Path to the folder with all the csv files for eg-footprint pixels.')
parser.add_option('--nside', dest='nside', type='int',
                  help='HEALPix resolution parameter. Default: 256',
                  default=256)
parser.add_option('--yr_cut', dest='yr_cut', type='int',
                  help='Year cut: 1, 3, 6, 10')
parser.add_option('--band', dest='band', type='str',
                  help='Filter band')
#######################################################################################
start_time = time.time()
# get the inputs
(options, args) = parser.parse_args()
print('## %s' % options)
outdir = options.outdir
db_path = options.db_path
eg_path = options.eg_path
nside = options.nside
yr_cut = options.yr_cut
band = options.band

# quick check
if yr_cut not in [1, 3, 6, 10]:
    raise ValueError('Currently can only work with yr_cut 1, 3, 6, 10. Inputs: %s'%yr_cut)

mag_cuts = {1: 24.75-0.1, 3: 25.35-0.1, 6: 25.72-0.1, 10: 26.0-0.1}

ilim = mag_cuts[yr_cut]

# set up the dbname
dbname = db_path.split('/')[-1].split('.db')[0]

# set up directory for bundle data
bundle_dir = '%s/bundle_data/' % outdir
os.makedirs(bundle_dir, exist_ok=True)
# set up directory for summary data
summary_dir = '%s/summary_data/' % outdir
os.makedirs(summary_dir, exist_ok=True)
# set up directory for log -- helps not re-running stuff.
log_dir = '%s/run_log/' % outdir
os.makedirs(log_dir, exist_ok=True)

# check if the db analysis has been run for this db already
log_filename = 'done_%s-band_y%s_ilim%s_nside%s.txt' % (band, yr_cut, ilim, nside)
if log_filename in os.listdir(log_dir):
    done_dbs = np.genfromtxt('%s/%s' % (log_dir, log_filename), dtype=str)
    if dbname in done_dbs:
        print('## Analysis already done for %s for y%s. Exiting.\n' % (dbname, yr_cut) )
        exit()
# -------------------------------------------------------------------------------------------------------------------------
opsdb = db.OpsimDatabase(db_path)

# set up the slicer
slicer = slicers.HealpixSlicer(lonCol='fieldRA', latCol='fieldDec',
                               latLonDeg=opsdb.raDecInDeg, nside=nside, useCache=False)
# set up the resultsdb object
resultsDb = db.ResultsDb(outDir=outdir)

# figure out the sql constraint
sqlconstraint = 'night <= %s and filter=="%s"'%(yr_cut * 365.25, band)
sqlconstraint += ' and note not like "DD%"'

# set up the dust map
dustmap = maps.DustMap(nside=nside, interp=False)

# set up the metric
metric = metrics.ExgalM5(lsstFilter=band)

# setup the bundle
bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, mapsList=[dustmap])
# set up the group.
grp = metricBundles.MetricBundleGroup({0: bundle}, opsdb, outDir=outdir,
                                      resultsDb=resultsDb, saveEarly=False)
grp.runAll()

# read in the footprint pixels
file = [f for f in os.listdir(eg_path) if f.__contains__(dbname) and \
                                            f.__contains__('_yr%s' % yr_cut)
       ][0]

if not file.__contains__('%s' % nside):
    raise ValueError('need eg-mask for nside %s' % (nside))

print('## reading in %s' % file)
eg_mask = np.array( np.genfromtxt('%s/%s' % (eg_path, file)) )
# find the indices for the in-survey pixels
good_pix = np.where( eg_mask == 0)[0]
# calculate statistics
area = len( good_pix ) * hp.nside2pixarea(nside=nside, degrees=True)
med_depth = np.median( bundle.metricValues.data[good_pix] )
std_depth = np.std( bundle.metricValues.data[good_pix] )
# set up what to write
to_write = '%s,%.2f,%.2f,%.2f\n' % ( dbname, area, med_depth, std_depth )
#
filename = 'median-depth_%s-band_eg-pixels_y%s_limi%s_nside%s.csv' % ( band, yr_cut, ilim, nside )
if filename not in os.listdir( summary_dir ):
    # file does not exist already so need to set up the header
    to_write = 'dbname,Area (deg2),%s-band depth: median,%s-band depth: std\n%s' % ( band, band, to_write )

# write to the file
txt_file = open('%s/%s'%(summary_dir, filename), 'a')
txt_file.write(to_write)
txt_file.close()

# save the depth map
outfile = 'exgal-depth-non-dd_%s_%s-band_nside%s.npz' % ( dbname, band, nside )
bundle.slicer.writeData('%s/%s'%(bundle_dir, outfile),
                        bundle.metricValues,
                        metricName = bundle.metric.name,
                        simDataName = bundle.runName,
                        constraint = bundle.constraint,
                        metadata = bundle.metadata,
                        displayDict = bundle.displayDict,
                        plotDict = bundle.plotDict
                       )

# add dbname to the log file
txt_file = open('%s/%s'%(log_dir, log_filename), 'a')
txt_file.write('%s \n' % dbname)
txt_file.close()
# all done.
print('## Time taken for %s: %.2f min\n' % (dbname, (time.time() - start_time )/60 ) )
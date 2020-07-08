##############################################################################################################
# This script runs the egFootprintMetric with the option to get the i-band depth map, calculates eg footprint
# area and depth statistics, and saves them. Also saves the depth data bundle.
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
import lsst.sims.maf.stackers as mafStackers

# new metric
from mafContrib.lssmetrics.depthLimitedNumGalMetric import depthLimitedNumGalMetric

#######################################################################################
from optparse import OptionParser
parser = OptionParser()
parser.add_option('--outdir', dest='outdir',
                  help='Path to the folder where all the output should be stored; should already exist.')
parser.add_option('--db_path', dest='db_path',
                  help='Path to the folder with the db to consider.')
parser.add_option('--nside', dest='nside', type='int',
                  help='HEALPix resolution parameter. Default: 256',
                  default=256)
parser.add_option('--yr_cut', dest='yr_cut', type='int',
                  help='Year cut: 1, 3, 6, 10')
parser.add_option('--wdith',
                  action='store_true', dest='wdith', default=False,
                  help='With translational dithers.')
parser.add_option('--older_cuts',
                  action='store_true', dest='older_cuts', default=False,
                  help='Use older depth cuts.')
parser.add_option('--fbs',
                  action='store_true', dest='fbs', default=False,
                  help='Treat as FBS output (no Proposal table).')
#######################################################################################
start_time = time.time()
# get the inputs
(options, args) = parser.parse_args()
print(options)
outdir = options.outdir
db_path = options.db_path
nside = options.nside
yr_cut = options.yr_cut
wdith = options.wdith
older_cuts = options.older_cuts
fbs = options.fbs
# quick check
if yr_cut not in [1, 3, 6, 10]:
    raise ValueError('Currently can only work with yr_cut 1, 3, 6, 10. Inputs: %s'%yr_cut)

# params for depth and ebv cuts
if older_cuts:
    mag_cuts = {1: 24.5, 3: 25.0, 6: 25.5, 10: 26.0}
else:
    #mag_cuts = {1: 24.75, 3: 25.35, 6: 25.72, 10: 26.0}
    mag_cuts = {1: 24.75-0.1, 3: 25.35-0.1, 6: 25.72-0.1, 10: 26.0-0.1}
lim_mag_i_ptsrc = mag_cuts[yr_cut]
nfilters_needed = 6
lim_ebv = 0.2

redshift_bin = '0.66<z<1.0'

# set up the dbname
dbname = db_path.split('/')[-1].split('.db')[0]

# set up the filter for which to get the depth data.
band = 'i'

# set up for dither vs. not.
if wdith:
    dither = 'RandomDitherPerNight'
else:
    dither = 'nodither'

# set up directory for bundle data
bundle_dir = '%s/bundle_data/' % outdir
os.makedirs(bundle_dir, exist_ok=True)
# set up directory for summary data
summary_dir = '%s/summary_data/' % outdir
os.makedirs(summary_dir, exist_ok=True)
# set up directory for log
log_dir = '%s/run_log/' % outdir
os.makedirs(log_dir, exist_ok=True)

# check if the db analysis has been run for this db already
log_filename = 'ngal_done_%s_y%s_lim%s%s_%s_nside%s.txt' % ( redshift_bin, yr_cut, band, lim_mag_i_ptsrc, dither.lower(), nside)
if log_filename in os.listdir(log_dir):
    done_dbs = np.genfromtxt('%s/%s' % (log_dir, log_filename), dtype=str)
    if dbname in done_dbs:
        raise ValueError('Analysis already done for %s for y%s.\n' % (dbname, yr_cut) )

# -------------------------------------------------------------------------------------------------------------------------
opsdb = db.OpsimDatabase(db_path)

# set up for dither vs. not.
if wdith:
    lonCol, latCol = 'randomDitherPerNightRa', 'randomDitherPerNightDec'
    stackerList = [mafStackers.RandomDitherPerNightStacker(degrees=opsdb.raDecInDeg, randomSeed=1000)]
else:
    lonCol, latCol = 'fieldRA', 'fieldDec'
    stackerList = []

# set up the slicer
slicer = slicers.HealpixSlicer(lonCol=lonCol, latCol=latCol,
                               latLonDeg=opsdb.raDecInDeg, nside=nside, useCache=False)
# set up the resultsdb object
resultsDb = db.ResultsDb(outDir=outdir)

# figure out the sql constraint
if fbs:
    sqlconstraint = 'night <= %s'%(yr_cut * 365.25)
    sqlconstraint += ' and note not like "DD%"'
else:
    propIds, propTags = opsdb.fetchPropInfo()
    wfdWhere = opsdb.createSQLWhere('WFD', propTags)
    sqlconstraint = '%s and night <= %s'%(wfdWhere, yr_cut * 365.25)

# set up the dust map
dustmap = maps.DustMap(nside=nside, interp=False)

# set up the metric
metric = depthLimitedNumGalMetric(nside=nside, filterBand=band,
                                  redshiftBin=redshift_bin, 
                                  nfilters_needed=nfilters_needed,
                                  lim_mag_i_ptsrc=lim_mag_i_ptsrc, lim_ebv=lim_ebv)
        
# setup the bundle
bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                    stackerList=stackerList, mapsList=[dustmap])
# set up the group.
grp = metricBundles.MetricBundleGroup({0: bundle}, opsdb, outDir=outdir,
                                      resultsDb=resultsDb, saveEarly=False)
grp.runAll()

# calculate the footprint area, etc.
good_pix = np.where(bundle.metricValues.mask == False)[0]
to_write = '%s,%.2e\n' % ( dbname, np.sum( bundle.metricValues.data[good_pix] ) )
#
filename = 'ngal_in_eg_stats_%s_y%s_lim%s%s_%s_nside%s.csv' % ( redshift_bin, yr_cut, band, lim_mag_i_ptsrc, dither.lower(), nside )
if filename not in os.listdir( summary_dir ):
    # file does not exist already so need to set up the header
    to_write = 'dbname,Ngal\n%s' % ( to_write )

# write to the file
txt_file = open('%s/%s'%(summary_dir, filename), 'a')
txt_file.write(to_write)
txt_file.close()

# save the depth map
outfile = 'ngal_in_eg_%s_%s_%s-band_lim%s%s_%s_nside%s.npz' % ( redshift_bin, dbname, band, band, lim_mag_i_ptsrc, dither.lower(), nside )
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
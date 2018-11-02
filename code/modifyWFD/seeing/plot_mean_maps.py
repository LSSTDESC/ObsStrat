# Goal here is to produce maps for specified quantities baseline2018a and
# pontus_2002 (wider footprint) after 1yr or 10-yrs.
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
########################################################################
from optparse import OptionParser
parser = OptionParser()
parser.add_option('--nside', dest='nside', type='int',
			help='HEALPix resolution parameter.', default=256)
parser.add_option('--one_yr',
                  action= 'store_true', dest= 'one_yr', default=False,
                  help= 'Run analyis for 1yr only.')
parser.add_option('--db_path', dest='db_path',
                  help='Path to OpSim db to consider.',
                  default='/global/cscratch1/sd/awan/dbs_wp_unzipped/')
parser.add_option('--outDir', dest='outDir',
                  help='Path to the folder where all the output should be stored. The folder must already exist.',
                  default='/global/cscratch1/sd/awan/lsst_output/seeing/')
parser.add_option('--specific_db', dest='specific_db',
                  help='Specific db to run the analysis for.',
                  default=None)
parser.add_option('--baseline_and_wide_only',
                  action='store_true', dest='baseline_and_wide_only', default=False,
                  help= 'baseline_and_wide_only')
parser.add_option('--quantity', dest='quantity',
                  help='Mean of the quantity to plot in each HEALPix pixel.: seeingFwhmEff, fiveSigmaDepth.',
                  default='seeingFwhmEff')

(options, args) = parser.parse_args()
one_yr = options.one_yr
dbs_path = options.db_path
outDir = options.outDir
specific_db = options.specific_db
baseline_and_wide_only = options.baseline_and_wide_only
nside = options.nside
quantity = options.quantity

########################################################################
dbfiles = [f for f in os.listdir(dbs_path) if f.endswith('db')]
if specific_db is not None:
    if specific_db in dbfiles:
        dbfiles = [specific_db]
    else:
        raise ValueError('%s is not in the dbs_path directory.'%specific_db)

if baseline_and_wide_only:
    dbfiles = [f for f in dbfiles if ((f.__contains__('baseline2018a') \
                                       or (f.__contains__('pontus_2002'))))]
dbfiles = sorted(dbfiles)
print('Running over %s'%dbfiles)

# set up
resultsDb = db.ResultsDb(outDir=outDir)

avgSeeingBundle = {}
bands = ['u', 'g', 'r', 'i', 'z', 'y']

########################################################################
# set things up for 1yr vs. 10yr
if one_yr:
    night = "night < 365.25"
    year_tag = "1yr"
else:
    night = None
    year_tag = "fullsurvey"
########################################################################
for dbfile in dbfiles:
    print('\n## Running %s for %s\n'%(dbfile, year_tag))
    startTime = time.time()
    # OpSim database
    opsdb = db.OpsimDatabase('%s/%s'%(dbs_path, dbfile))
    # WFD only
    propIds, propTags = opsdb.fetchPropInfo()
    wfdWhere = opsdb.createSQLWhere('WFD', propTags)

    # implement per night dithers
    stackerList = [mafStackers.RandomDitherPerNightStacker(degrees=opsdb.raDecInDeg, randomSeed=1000)]
    slicer= slicers.HealpixSlicer(lonCol='randomDitherPerNightRa', latCol='randomDitherPerNightDec',
                                  latLonDeg=opsdb.raDecInDeg, nside=nside, useCache=False)
    meanMetric= metrics.MeanMetric(col=quantity)   # for avg quantity per HEALpix pixel

    db_key = dbfile.split('.db')[0]
    avgSeeingBundle[db_key] = {}

    for band in bands:
        if wfdWhere is None and night is None:
            sqlconstraint = 'filter=="%s"'%(band)
        elif wfdWhere is not None and night is None:
            sqlconstraint = '%s and filter=="%s"'%(wfdWhere, band)
        elif wfdWhere is None and night is not None:
            sqlconstraint = '%s and filter=="%s"'%(night, band)
        else:
            sqlconstraint = '%s and %s and filter=="%s"'%(wfdWhere, night, band)
        avgSeeingBundle[db_key][band] = metricBundles.MetricBundle(meanMetric, slicer, sqlconstraint,
                                                                   stackerList= stackerList)

    sGroup = metricBundles.MetricBundleGroup(avgSeeingBundle[db_key], opsdb, outDir=outDir,
                                             resultsDb=resultsDb,saveEarly= False)
    sGroup.runAll()
    print('\n## %s: Time taken: %.2f min\n'%(db_key, (time.time()-startTime)/60.))

# plot skymaps
for band in bands:
        plt.clf()
        fig, axes= plt.subplots(1,2)
        db_tag = ''
        for i, db_key in enumerate(avgSeeingBundle.keys()):
            db_tag += '_%s'%db_key

            if (i==0): # define the color range
                inSurveyIndex = np.where(avgSeeingBundle[db_key][band].metricValues.mask == False)[0]
                median = np.median(avgSeeingBundle[db_key][band].metricValues.data[inSurveyIndex])
                stddev = np.std(avgSeeingBundle[db_key][band].metricValues.data[inSurveyIndex])

                nTicks = 5
                if quantity=='seeingFwhmEff':
                    colorMin = 0.8
                    colorMax = 1.4
                elif quantity=='fiveSigmaDepth':
                    if band=='u':
                        colorMin, colorMax = 22.7, 23.4
                    elif band=='g':
                        colorMin, colorMax = 24.1, 24.75
                    elif band=='r':
                        colorMin, colorMax = 23.5, 24.5
                    elif band=='i':
                        colorMin, colorMax = 23.25, 23.75
                    elif band=='z':
                        colorMin, colorMax = 22.4, 23.0
                    elif band=='y':
                        colorMin, colorMax = 21.5, 22.2
                    else:
                        raise ValueError('Somethings wrong.')
                else:
                    colorMin = median-2.5*stddev
                    colorMax = median+2.5*stddev

                increment = (colorMax-colorMin)/float(nTicks)
                ticks = np.arange(colorMin+increment, colorMax, increment)

            plt.axes(axes[i])
            hp.mollview(avgSeeingBundle[db_key][band].metricValues.filled(avgSeeingBundle[db_key][band].slicer.badval),
                        flip='astro', rot=(0,0,0) , hold=True,
                        min=colorMin, max=colorMax,
                        title= '%s: WFD'%db_key, cbar=False)
            hp.graticule(dpar=20, dmer=20, verbose=False)

        im = plt.gca().get_images()[0]
        cbaxes = fig.add_axes([0.25, 0.38, 0.5, 0.01]) # [left, bottom, width, height]
        cb = plt.colorbar(im, orientation='horizontal',
                          ticks=ticks, format='%.2f', cax=cbaxes)
        cb.set_label('%s-band: %s after %s'%(band, quantity, year_tag), fontsize=14)
        cb.ax.tick_params(labelsize=14)

        fig.set_size_inches(15,15)
        filename = '%s%s_%s_nside%s_%sband.png'%(quantity, db_tag, year_tag, nside, band)
        plt.savefig('%s/%s'%(outDir, filename), format='png',  bbox_inches='tight')
        print('## Saved %s in outDir.\n'%filename)
        #plt.show()

# plot histograms
colors = ['g', 'r', 'k', 'b', 'm', 'c']
plt.clf()
for i, band in enumerate(bands):
    db_tag = ''
    for db_key in avgSeeingBundle.keys():
        db_tag += '_%s'%db_key

        if quantity=='seeingFwhmEff':
            bins = np.arange(0.5, 2, 0.01)
        elif quantity=='fiveSigmaDepth':
            bins = np.arange(20, 25, 0.01)
        else:
            bins = None

        # set the linestyle
        if db_key.__contains__('baseline'):
            sty = 'solid'
        else:
            sty = 'dashed'

        plt.hist(avgSeeingBundle[db_key][band].metricValues.filled(avgSeeingBundle[db_key][band].slicer.badval),
                 bins=bins, label='%s-band: %s'%(band, db_key),
                 histtype='step', color=colors[i], linestyle=sty)
plt.ylabel('Counts')
plt.xlabel('%s %s'%(year_tag, quantity))
plt.legend(bbox_to_anchor=(1,1))
plt.gcf().set_size_inches(10,5)
filename = '%s%s_histogram_%s_nside%s.png'%(quantity, db_tag, year_tag, nside)
plt.savefig('%s/%s'%(outDir, filename), format='png',  bbox_inches='tight')
print('## Saved %s in outDir.\n'%filename)
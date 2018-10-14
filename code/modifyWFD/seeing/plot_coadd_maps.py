# Goal here is to read-in the already saved data for coadded depth maps for kraken_2026 (baseline)
# and pontus_2002 (wider footprint) after 1yr and 10-yrs and plot the skymaps side-by-side.
#
# Humna Awan: humna.awan@rutgers.edu
#
########################################################################
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
import os
import healpy as hp
import copy
from collections import OrderedDict
import lsst.sims.maf.metricBundles as metricBundles   # need MAF installed; for reading in coadd data

########################################################################
from optparse import OptionParser
parser = OptionParser()
parser.add_option('--data_dir', dest='data_dir',
                  help='Path to directory with all the coadd outputs.',
                  default='/global/cscratch1/sd/awan/lsst_output/coadd_output_allwps_perNight')
parser.add_option('--outDir', dest='outDir',
                  help='Path to the folder where all the output should be stored. The folder must already exist.',
                  default='/global/cscratch1/sd/awan/lsst_output/coadd_maps')
parser.add_option('--one_yr',
                  action= 'store_true', dest= 'one_yr', default=False,
                  help= 'Run analyis for 1yr only.')
parser.add_option('--noMWdust',
                  action='store_true', dest='noMWdust', default=False,
                  help= 'Use the tag to not include MW dust extinction.')
parser.add_option('--owsee_tag', dest='owsee_tag',
                  help='Tag to identify an owsee outputs',
                  default=None)

(options, args) = parser.parse_args()
print(options)
one_yr = options.one_yr
includeDustExtinction = not options.noMWdust
data_dir = options.data_dir
outDir = options.outDir
owsee_tag = options.owsee_tag

########################################################################################################################
# set some things up
#data_dir = '/global/homes/a/awan/LSST/output/coadd_output'   # where coadd depth maps are
#outDir = '/global/homes/a/awan/LSST/output/coadd_maps'

orderBand = ['u', 'g', 'r', 'i', 'z', 'y']
nside = 256

print('data_dir: %s'%data_dir)
print('outDir: %s'%outDir)
print('nside: %s\n'%nside)

if one_yr:
    file_yearTag = '1yearCut'
    plt_yearTag = '1yr'
else:
    file_yearTag = 'fullSurveyPeriod'
    plt_yearTag = 'fullsurvey'

if includeDustExtinction:
    dust_tag = 'withDustExtinction'
else:
    dust_tag = 'noDustExtinction'

dbs = os.listdir(data_dir)
dbs = [f for f in dbs if f.__contains__('zBand') and f.__contains__(file_yearTag) and \
               (f.__contains__('baseline2018a') or f.__contains__('pontus_2002'))]

dbs = [f.split('zBand_')[-1].split('_%s'%file_yearTag)[0] for f in dbs]
if owsee_tag is not None:
    dbs = [f for f in dbs if f.__contains__('%s_'%owsee_tag)]
dbs = sorted(dbs)
print('Creating plots for %s'%dbs)

########################################################################################################################
data_bundle = OrderedDict()
for band in orderBand:
    for db in dbs:
        folder = 'coaddM5Analysis_nside%s_%s_0pixelRadiusForMasking_%sBand_%s_%s_directory/'%(nside, dust_tag, band, db,
                                                                                              file_yearTag)

        path = '%s/%s/unmaskedCoaddData/'%(data_dir, folder)

        filenames = [f for f in os.listdir(path) if f.endswith('.npz')]
        print('Reading %s from\n%s/unmaskedCoaddData.\n'%(filenames, folder))

        for filename in filenames:
            mB = metricBundles.createEmptyMetricBundle()
            mB.read('%s/%s'%(path, filename))
            data_bundle['%s_%s'%(db, band)]= mB

# plot skymaps
for band in orderBand:
    plt.clf()
    fig, axes= plt.subplots(1,2)
    db_tag = ''     # for the filename
    for i, db_key in enumerate(dbs):
        key = '%s_%s'%(db_key, band)
        db_tag += '_%s'%db_key
        if i==0:
            #inSurveyIndex = np.where(data_bundle[key].metricValues.mask == False)[0]
            #median = np.median(data_bundle[key].metricValues.data[inSurveyIndex])
            #stddev = np.std(data_bundle[key].metricValues.data[inSurveyIndex])
            nTicks = 5
            #colorMin = median-2.5*stddev
            #colorMax = median+2.5*stddev
            if band=='u':
                colorMin, colorMax = 23.5, 25.0
            elif band=='g':
                colorMin, colorMax = 24.5, 26.5
            elif band=='r':
                colorMin, colorMax = 25.0, 26.5
            elif band=='i':
                colorMin, colorMax = 24.25, 26.25
            elif band=='z':
                colorMin, colorMax = 23.5, 25.25
            elif band=='y':
                colorMin, colorMax = 23.0, 24.25
            else:
                raise ValueError('Somethings wrong.')
            if not one_yr:
                colorMin += 1.0
                colorMax += 1.0
            increment = (colorMax-colorMin)/float(nTicks)
            ticks = np.arange(colorMin+increment, colorMax, increment)
        plt.axes(axes[i])
        hp.mollview(data_bundle[key].metricValues.filled(data_bundle[key].slicer.badval),
                    flip='astro', rot=(0,0,0) , hold= True,
                    min=colorMin, max=colorMax,
                    title= '%s: WFD'%db_key, cbar=False)
        hp.graticule(dpar=20, dmer=20, verbose=False)
    im = plt.gca().get_images()[0]
    cbaxes = fig.add_axes([0.25, 0.38, 0.5, 0.01]) # [left, bottom, width, height]
    cb = plt.colorbar(im,  orientation='horizontal',
                      ticks=ticks, format='%.2f', cax=cbaxes)
    cb.set_label(r'%s-band: 5$\sigma$ Coadded Depth after %s; %s'%(band, plt_yearTag, dust_tag), fontsize=14)
    cb.ax.tick_params(labelsize=14)
    fig.set_size_inches(15,15)
    filename = 'fiveSigCoaddDepth_%s%s_%s_nside%s_%sband.png'%(dust_tag, db_tag, plt_yearTag, nside, band)
    plt.savefig('%s/%s'%(outDir, filename), format='png',  bbox_inches='tight')
    print('## Saved %s in outDir.\n'%filename)
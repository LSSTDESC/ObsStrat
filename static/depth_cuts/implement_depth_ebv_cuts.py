# Goal here is to implement various depth (+ an extinction cut) to define the footprint for extragalactic science.
#
# Need the coadded depth data for various year cuts (e.g., 1yr, 3yr, 10yr); calculated by bash_calc_coadd_all.sh
#
# Humna Awan: humna.awan@rutgers.edu
#
########################################################################################################################
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import healpy as hp
import copy
from collections import OrderedDict
import lsst.sims.maf.metricBundles as metricBundles   # need MAF installed; for reading in coadd data
# for galactic latitude histogram
from astropy.coordinates import SkyCoord
from astropy import units as u
# for EBV map from MAF
import lsst.sims.maf.db as db
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.stackers as mafStackers   # stackers in sims_maf
import lsst.sims.maf.maps as maps
import time
import pickle
########################################################################################################################
from optparse import OptionParser
parser = OptionParser()
parser.add_option('--nside', dest='nside', type='int',
			help='HEALPix resolution parameter.', default=256)
parser.add_option('--db_path', dest='db_path',
                  help='Path to OpSim db to consider.',
                  default='/global/cscratch1/sd/awan/dbs_wp_unzipped/baseline2018a.db')
parser.add_option('--coadd_data_dir', dest='coadd_data_dir',
                  help='Path to the folder where coadd depth maps are.',
                  default='/global/homes/a/awan/LSST/output/coadd_output_allwps_perNight/')
parser.add_option('--mag_cuts', dest='mag_cuts',
                  help='List of depth cuts to consider',
                  default='22.0, 23.0, 24.0, 24.5, 25.0, 25.3, 25.5, 26.0, 26.5')
parser.add_option('--yr_cuts', dest='yr_cuts',
                  help='Years cuts to consider, e.g. 1yr, 10yr.',
                  default= '1yr, 10yr')
parser.add_option('--chosen_cuts', dest='chosen_cuts',
                  help='Finalized cuts for the yr_cuts, in the same as order as yr_cuts input.',
                  default= '24.5, 26.0')
parser.add_option('--ebv_cut',
                  action='store_true', dest='ebv_cut', default=False,
                  help= 'Set to True if want to implement a EBV cut (discard pixels with EBV>0.2)')
parser.add_option('--save_stuff',
                  action='store_true', dest='save_stuff', default=False,
                  help= 'Set to True if want to save the plots and the readme etc.')
parser.add_option('--dont_show_plots',
                  action='store_true', dest='dont_show_plots', default=False,
                  help= 'Set to True if dont want to show the plots.')
parser.add_option('--outDir', dest='outDir',
                  help='Path to the folder where all the output should be stored; directory should exist already.',
                  default='/global/homes/a/awan/LSST/output/')
parser.add_option('--outDir_md', dest='outDir_md',
                  help='Path to the folder where md files should be stored; directory should exist already. If None, same as outDir.',
                  default=None)
parser.add_option('--debug',
                  action='store_true', dest='debug', default=False,
                  help= 'Set to True if want to debug: basically run the analysis as in DESC SRD v1.')
parser.add_option('--slair',
                  action='store_true', dest='slair', default=False,
                  help= 'Use the tag to specify a slair run.')
parser.add_option('--bands', dest='bands',
                  help='Bands to consider.',
                  default= 'u, g, r, i, z, y')
########################################################################################################################
startTime = time.time()
(options, args) = parser.parse_args()
print('\nOptions: %s'%options)

# read in the inputs
nside = options.nside
dbpath = options.db_path
data_dir = options.coadd_data_dir
mag_cuts = options.mag_cuts
yr_cuts = options.yr_cuts
final_mag_cuts = options.chosen_cuts
ebv_cut = options.ebv_cut
save_stuff = options.save_stuff
dont_show_plots = options.dont_show_plots
outDir = options.outDir
outDir_md = options.outDir_md
slair = options.slair
bands = options.bands

# format the mag cuts
mag_cuts = [float(f) for f in list(mag_cuts.split(','))]
# format the yr cuts
yr_cuts = [f.strip() for f in list(yr_cuts.split(','))]
print('\nYear cuts to be implemented: %s'%(yr_cuts))
# format the final mag cuts
cuts = [float(f) for f in list(final_mag_cuts.split(','))]
if len(cuts)!= len(yr_cuts):
    raise ValueError('Need a final cut for each of yr_cuts. Currently have %s for %s.'%(cuts, yr_cuts))
# ensure that the final cut is included in the mag_cut list. if not, add it.
for cut in cuts:
    if cut not in mag_cuts:
        mag_cuts.append(cut)
        print('Adding %s cut to the mag_cut list since it is one of the final cuts.'%cut)
chosen_cuts = {}
for i, yr_cut in enumerate(yr_cuts):
    chosen_cuts[yr_cut] = cuts[i]
print('\nFinal depth cuts to be implemented: %s'%(chosen_cuts))
# figure out the outDir for md files.
if outDir_md is None: outDir_md = outDir
# set up for debugging
debug = options.debug
if debug:
    dbpath = '/global/cscratch1/sd/awan/dbs_old_unzipped/minion_1016_sqlite_new_dithers.db'
    data_dir = '/global/homes/a/awan/LSST/output/coadd_output_minion1016/'

# figure out the band
if bands.__contains__(','):
    orderBand = list(bands.split(','))
else:
    orderBand = bands
########################################################################
# set some things up
dbname = dbpath.split('/')[-1].split('.db')[0]

print('\ndata_dir: %s'%data_dir)
print('dbname: %s'%dbname)
print('nside: %s\n'%nside)

########################################################################################################################
# get the files and put the data in a bundle
print('## Reading in the data ... \n')
data_bundle = OrderedDict()
for yr_cut in yr_cuts:
    if yr_cut=='10yr': file_yearTag = 'fullSurveyPeriod'
    else: file_yearTag = '%syearCut'%(yr_cut.split('yr')[0])
        
    for band in orderBand:
        folder = 'coaddM5Analysis_nside%s_withDustExtinction_0pixelRadiusForMasking_%sBand_%s_%s_directory/'%(nside, band,
                                                                                                              dbname,
                                                                                                              file_yearTag)
        path = '%s/%s/unmaskedCoaddData/'%(data_dir, folder)
        filenames = [f for f in os.listdir(path) if f.endswith('.npz')]
        print('Reading %s from\n%s/unmaskedCoaddData.\n'%(filenames, folder))
        
        if len(filenames)>1:
            raise ValueError('Have more than one npz file for %s band for %s data: %s'%(band, yr_cut, filenames))
        else:
            dither = filenames[0].split('%s_'%band)[-1]
            dither = dither.split('.npz')[0]
            mB = metricBundles.createEmptyMetricBundle()
            mB.read('%s/%s'%(path, filenames[0]))
            data_bundle['%s_%s'%(yr_cut, band)]= mB

print('\n%s dithers'%dither)
########################################################################################################################
# check the improvement factor between 1, 10yr data
if yr_cuts.__contains__('1yr') and yr_cuts.__contains__('10yr'):
    print('\n## Calculating improvement in fluxes between 1yr and 10yrs ... ')
    allImprovs = []
    for band in orderBand:
        yr_cut = '1yr_%s'%band
        in_survey_positive = np.where((data_bundle[yr_cut].metricValues.mask == False) & \
                                     (data_bundle[yr_cut].metricValues.data > 0))[0]
        one_yr_med = np.median(data_bundle[yr_cut].metricValues.data[in_survey_positive])

        yr_cut = '10yr_%s'%band
        in_survey_positive = np.where((data_bundle[yr_cut].metricValues.mask == False) & \
                                     (data_bundle[yr_cut].metricValues.data > 0))[0]
        ten_yr_med = np.median(data_bundle[yr_cut].metricValues.data[in_survey_positive])

        one_yr_flux = 10**(-one_yr_med/2.5)
        ten_yr_flux = 10**(-ten_yr_med/2.5)
        print('%s-band: improvement factor in flux: %s'%(band, one_yr_flux/ten_yr_flux))
        allImprovs.append(one_yr_flux/ten_yr_flux)

    print('Wanted improvement factor over ten years: %s'%np.sqrt(10.))
    print('Mean improvement factor across ugrizy: %s'%np.mean(allImprovs))

##
if debug:
    if not (yr_cuts.__cotains__('1yr') and yr_cuts.__contains__('10yr')):
        raise ValueError('Need 1yr, 10yr data to run the debug analysis.')
    # The improvement from 1-10yr is too good -- 1yr in minion1016 strongly prefers DDFs so 1yr defined by the number
    # of the nights doesnt have 10% of total visits (has <7%). For now, renormalize the 1yr depth s.t. median matches
    # that after 10% of 10-year WFD observations.
    # As implemented in https://github.com/humnaawan/MAF-Related-Notebooks/blob/master/DESC-SRD/srd_depth_cuts.py
    wanted1yr_medianDepth= {'g': 25.377165833786055, 'i': 24.910057884620223, 'r': 25.565945074516804,
                            'u': 23.795160853950424, 'y': 23.315667199085482, 'z': 24.002597276614527}

    print('\n## Renormalizing 1yr depth data ... ')
    inSurveyIndex = {}
    for key in data_bundle.keys():
        inSurveyIndex[key] = np.where(data_bundle[key].metricValues.mask == False)[0]
        if key.__contains__('1yr'):
            band = key.split('1yr_')[1]
            print(band)
            medDepth = np.median(data_bundle[key].metricValues.data[inSurveyIndex[key]])
            print('Median depth as read: %s'%np.median(data_bundle[key].metricValues.data[inSurveyIndex[key]]))
            delm = wanted1yr_medianDepth[band]-medDepth
            print('m_wanted-m_current = %s'%delm)
            data_bundle[key].metricValues.data[:] += delm
            print('Renormalized map. \nNew median: %s\n'%np.median(data_bundle[key].metricValues.data[inSurveyIndex[key]]))

    # re-check the improvement factor between 1, 10yr data
    print('## Re-calculating improvement in fluxes between 1yr and 10yrs after the renormalizing ... ')
    allImprovs = []
    for band in orderBand:
        yr_cut = '1yr_%s'%band
        in_survey_positive = np.where((data_bundle[yr_cut].metricValues.mask == False) & \
                                      (data_bundle[yr_cut].metricValues.data > 0))[0]
        one_yr_med = np.median(data_bundle[yr_cut].metricValues.data[in_survey_positive])

        yr_cut = '10yr_%s'%band
        in_survey_positive = np.where((data_bundle[yr_cut].metricValues.mask == False) & \
                                      (data_bundle[yr_cut].metricValues.data > 0))[0]
        ten_yr_med = np.median(data_bundle[yr_cut].metricValues.data[in_survey_positive])

        one_yr_flux = 10**(-one_yr_med/2.5)
        ten_yr_flux = 10**(-ten_yr_med/2.5)
        print('%s-band: improvement factor in flux: %s'%(band, one_yr_flux/ten_yr_flux))
        allImprovs.append(one_yr_flux/ten_yr_flux)

    print('\nWanted improvement factor over ten years: %s'%np.sqrt(10.))
    print('Mean improvement factor across ugrizy: %s'%np.mean(allImprovs))
else:
    inSurveyIndex = {}
    for key in data_bundle.keys():
        inSurveyIndex[key] = np.where(data_bundle[key].metricValues.mask == False)[0]

########################################################################################################################
########################################################################################################################
# calculate some stats
areaPerPixel= hp.pixelfunc.nside2pixarea(nside=nside, degrees=True)

def calc_stats(bundle, index, allBandInds=False, return_stuff=False):
    # index must have the same keys as bundle
    if (bundle.keys()!=index.keys()) and not allBandInds:
        raise ValueError('index must have the same keys as bundle:\n%s\n%s'%(bundle.keys(), index.keys()))
        
    if return_stuff and allBandInds: 
        stuff_to_return = {}
        for key in ['5$\sigma$ Depth: Median', '5$\sigma$ Depth: Std', 'Area (deg$^2$)']:
            stuff_to_return[key] = {}

    md_print = ''
    header, sep = '| - ', '| ----:'
    med_depth, std_depth, area = '| 5$\sigma$ Depth: Median ', '| 5$\sigma$ Depth: Std ', '| Area (deg$^2$) '
    yr = None
    for key in bundle:
        if yr is None: yr = key.split('yr')[0]+'yr'
        
        current_yr = key.split('yr')[0]+'yr'
        if current_yr!=yr:
            md_print += '%s\n%s\n%s\n%s\n%s\n\n'%(header, sep, med_depth, std_depth, area)
            header, sep = '| - ', '| ----:'
            med_depth, std_depth, area = '| 5$\sigma$ Depth: Median ', '| 5$\sigma$ Depth: Std ', '| Area (deg$^2$) '
            yr = current_yr
        
        if allBandInds: index_key = current_yr
        else: index_key = key
            
        med = np.nanmedian(bundle[key].metricValues.data[index[index_key]])
        std = np.nanstd(bundle[key].metricValues.data[index[index_key]])
        sarea = (len(index[index_key])*areaPerPixel)
        
        if return_stuff and allBandInds:
            stuff_to_return['5$\sigma$ Depth: Median'][key] = med
            stuff_to_return['5$\sigma$ Depth: Std'][key] = std
            stuff_to_return['Area (deg$^2$)'][index_key] = sarea
            
        header += '| %s '%key
        sep += '|:----:'
        med_depth += '| %.2f '%med
        std_depth += '| %.2f '%std
        area += '| %.2f '%sarea
    
    md_print += '%s\n%s\n%s\n%s\n%s\n'%(header, sep, med_depth, std_depth, area)
    print(md_print)

    if return_stuff: return stuff_to_return, md_print
########################################################################################################################
########################################################################################################################
# Calculate stats in the survey region (unmasked; no constraints on depth, i.e., even have negative depths rn).
print('\n#### Stats: no constraints on depth, i.e., even have negative depths')
calc_stats(bundle=data_bundle, index=inSurveyIndex, allBandInds=False)

########################################################################################################################
# Find the area common to all-6 bands with depths>0 in all.
allBandPixels = {}  # dictionary for pixels that are common in all six bands with depth>0.

for key in data_bundle:
    index = np.where((data_bundle[key].metricValues.mask == False) & \
                     (data_bundle[key].metricValues.data > 0))[0]
    # save the indices
    yrTag = key.split('yr')[0]+'yr'
    if yrTag not in allBandPixels.keys():
        allBandPixels[yrTag]= index  
    else:
        allBandPixels[yrTag]= list(set(allBandPixels[yrTag]).intersection(index))

for key in allBandPixels: allBandPixels[key] = np.array(allBandPixels[key])
########################################################################################################################
# Calculate the stats for the all-band footprint
print('\n#### Stats: considering area common to all-6 bands with depths>0 in all')
calc_stats(bundle=data_bundle, index=allBandPixels, allBandInds=True)

########################################################################################################################
# implement depth cuts in i-band and save the pixel numbers
iCutPixels = {}
# run over different cuts
for mag_cut in mag_cuts:
    if mag_cut not in iCutPixels.keys():
        iCutPixels[mag_cut] = {}
        
    for yrTag in allBandPixels:
        if yrTag not in iCutPixels[mag_cut].keys():
            iCutPixels[mag_cut][yrTag] = {}
        
        # find the pixels satisfying the iBand cut.
        iBandCutInd = np.where((data_bundle['%s_i'%yrTag].metricValues.data[allBandPixels[yrTag]]>mag_cut))[0]
        iCutPixels[mag_cut][yrTag] = np.array(allBandPixels[yrTag])[iBandCutInd] # store

########################################################################################################################
# Calculate the stats in the survey region (unmasked; no constraints on depth, i.e., even have negative depths rn).
dat_keys = ['Area (deg$^2$)', '5$\sigma$ Depth: Median', '5$\sigma$ Depth: Std']

########################################################################################################################
# plots for area and depth variations as a funtion of mag cuts
# need to create a list of all the stats for cleaner plots.
stats_allmags = {}
for mag_cut in mag_cuts:
    print('\n#### Stats: i>%s in area common to all six bands with depths>0 in all'%mag_cut)
    stats, _ = calc_stats(bundle=data_bundle, index=iCutPixels[mag_cut],
                       allBandInds=True, return_stuff=True) # area for yr_cut; depth stuff for yrcut_band
    for dat_key in dat_keys:
        if dat_key not in stats_allmags: stats_allmags[dat_key] = {}
        for key in stats[dat_key].keys():
            if key.__contains__('_'):  # need to separate yrcut from yrcut_<band>
                sp = key.split('_')
                yr, band = sp[0], sp[1]
            else:
                yr, band = key, None
                
            if band is None: # all-band keys: yrcut
                if yr not in stats_allmags[dat_key]:
                    stats_allmags[dat_key][yr] = []
            else: # need to account for the band for each of the yr_cut
                if yr not in stats_allmags[dat_key]:
                    stats_allmags[dat_key][yr] = {}
                if band not in  stats_allmags[dat_key][yr]:
                     stats_allmags[dat_key][yr][band] = []
            
            if band is not None:
                stats_allmags[dat_key][yr][band].append(stats[dat_key][key])
            else:
                stats_allmags[dat_key][yr].append(stats[dat_key][yr])

print('\n## Plotting area and depth variations for mag_cuts: %s ...'%mag_cuts)
colors = ['m', 'g', 'b', 'r', 'k', 'c']
fontsize = 14
# plot
plt.clf()
nrows, ncols = len(yr_cuts), 3
fig, axes = plt.subplots(nrows, ncols)
fig.subplots_adjust(wspace=.2, hspace=.3)

for i, dat_key in enumerate(dat_keys):
    if dat_key.__contains__('Depth'):  # band-specific statistic
        for j, band in enumerate(stats_allmags[dat_key][yr_cuts[0]].keys()):
            for m, yr_cut in enumerate(yr_cuts):
                # plot yr_cut data for this statistic
                axes[m, i].plot(mag_cuts, stats_allmags[dat_key][yr_cut][band], 'o-',
                                color=colors[j], label='%s-band'%band)
    else:
        for m, yr_cut in enumerate(yr_cuts):
            # plot yr_cut area
            axes[m, i].plot(mag_cuts, stats_allmags[dat_key][yr_cut], 'o-')
for row in range(nrows):
    axes[row, 0].ticklabel_format(style='sci',scilimits=(-3,4), axis='y')  # area
    axes[row, 2].legend(loc='upper right', ncol=2, fontsize=fontsize-2)
    for col in range(ncols):
        axes[row, col].set_ylabel(dat_keys[col], fontsize=fontsize)    
        axes[row, col].tick_params(axis='x', labelsize=fontsize-2)
        axes[row, col].tick_params(axis='y', labelsize=fontsize-2)
# set all the ylims
ymin_final, ymax_final = [10**9, 10**9, 10**9], [0, 0, 0]
for i, col in enumerate(range(ncols)):
    for row in range(nrows):
        ymin, ymax = axes[row, col].get_ylim()
        ymin_final[i] = min(ymin_final[i], ymin)
        ymax_final[i] = max(ymax_final[i], ymax)
for i, col in enumerate(range(ncols)):
    for row in range(nrows):
        axes[row, col].set_ylim(ymin_final[i], ymax_final[i])
# set the title
for m, yr_cut in enumerate(yr_cuts):
    axes[m,1].set_title(yr_cut, fontsize=fontsize)

axes[m, 1].set_xlabel('i-band cut (i>?) (in all-band footprint with all depth > 0)', fontsize=fontsize)
fig.set_size_inches(20, nrows*5)
if save_stuff:
    filename = 'stats_variation_%s_nside%s_%s.png'%(dbname, nside, dither)
    plt.savefig('%s/%s'%(outDir, filename), format= 'png', bbox_inches='tight')
    print('## Saved %s in %s.'%(filename, outDir))
if not dont_show_plots:
    plt.show()
else:
    plt.close('all')

########################################################################################################################
# plot galactic latitude and EBV histograms for different cuts
print('\n## Plotting galactic latitude and EBV histograms for mag_cuts: %s'%mag_cuts)
# import EBV map from MAF
# read in the database
if slair:
    # slair database
    opsdb = db.Database(dbpath, defaultTable='observations')
    cols = ['RA', 'dec', 'night']
    raDecInDeg = True
else:
    # OpSim database
    opsdb = db.OpsimDatabase(dbpath)
    cols = ['fieldRA', 'fieldDec', 'night']
    raDecInDeg = opsdb.raDecInDeg
# decide on the columns to get.
#if debug:  # for minion1016
#    cols = ['fieldID', 'fieldRA', 'fieldDec', 'night']
#else:
#    cols = ['fieldId', 'fieldRA', 'fieldDec', 'night']

simdata = opsdb.fetchMetricData(cols, sqlconstraint=None)

# decide on the pointing RA, Dec
if slair:
    lonCol, latCol = 'RA', 'dec'
else:
    if dither.__contains__('NoDither'):  # for slair, altsched
        lonCol, latCol = 'fieldRA', 'fieldDec'
    else:
        dither_timescale = dither.split('Dither')[-1]
        lonCol, latCol = 'randomDither%sRa'%dither_timescale, 'randomDither%sDec'%dither_timescale
        # set up stacker
        if dither.__contains__('FieldPerNight'):
            s = mafStackers.RandomDitherFieldPerNightStacker(degrees=raDecInDeg, randomSeed=1000)
        elif dither.__contains__('PerNight'):
            s = mafStackers.RandomDitherPerNightStacker(degrees=raDecInDeg, randomSeed=1000)
        elif dither.__contains__('FieldPerVisit'):
            s = mafStackers.RandomDitherFieldPerVisitStacker(degrees=raDecInDeg, randomSeed=1000)
        else:
            raise ValueError('Unsure of what stacker to consdier for %s dithers.'%dither)
        simdata = s.run(simdata)

# slice the data
dustmap = maps.DustMap(nside=nside)
slicer = slicers.HealpixSlicer(lonCol=lonCol, latCol=latCol,
                               latLonDeg=raDecInDeg, nside=nside)
slicer.setupSlicer(simdata)
result = dustmap.run(slicer.slicePoints)
ebv_map = result['ebv']

########################################################################################################################
# histograms: latitude, extinction
bins_b = np.arange(-90, 90, 0.5)
bins_ebv = np.arange(-3.5, 1.0, 0.05)
colors = ['m', 'g', 'b', 'r', 'c', 'y']

# plot
plt.clf()
nrows, ncols = len(yr_cuts), 2
fig, axes = plt.subplots(nrows, ncols)
fig.subplots_adjust(wspace=0.2, hspace=0.3)

max_counts = 0  # for EBV histogram; needed for plotting constant EBV lines
for i, yr in enumerate(yr_cuts):
    linestyle = 'solid'
    # plot the galactic latitude histogram for no-cut
    lon, lat = hp.pix2ang(nside=nside, ipix=allBandPixels[yr], lonlat=True)
    c = SkyCoord(ra= lon*u.degree, dec=lat*u.degree)
    axes[i, 0].hist(c.galactic.b.degree, label='all-band; all depths>0',
                    bins=bins_b, histtype='step', lw=2, color='k')
    
    # plot the EBV histogram for no cut
    cts, _, _ = axes[i, 1].hist(np.log10(ebv_map[allBandPixels[yr]]), # label='all-band; all depths>0',
                                bins=bins_ebv, histtype='step', lw=2, color='k')
    max_counts = max(max_counts, max(cts))
    
    # now loop over mag cuts
    for j, mag_cut in enumerate(iCutPixels):
        if j>=len(colors): linestyle = 'dashed'  # out of colors so change linestyle
            
        # plot the galactic latitude histogram
        lon, lat= hp.pix2ang(nside=nside, ipix=iCutPixels[mag_cut][yr], lonlat= True)
        c = SkyCoord(ra= lon*u.degree, dec=lat*u.degree)
        axes[i, 0].hist(c.galactic.b.degree, label= 'i>%s'%mag_cut,
                        bins=bins_b, histtype='step', lw=2,
                        color=colors[j%len(colors)], linestyle=linestyle)
        
        # plot the EBV histogram
        cts, _, _ = axes[i, 1].hist(np.log10(ebv_map[iCutPixels[mag_cut][yr]]), #label= 'i>%s'%mag_cut,
                                    bins=bins_ebv, histtype='step', lw=2,
                                    color=colors[j%len(colors)], linestyle=linestyle)
        max_counts = max(max_counts, max(cts))

for row in range(nrows):
    x = np.arange(0, max_counts, 10)
    for ebv in [0.2, 0.3]:
        axes[row, 1].plot(np.log10([ebv]*len(x)), x, '-.', label='EBV: %s'%ebv)
    axes[row, 0].legend(bbox_to_anchor=(2.65, 1.0), fontsize=fontsize-2)
    axes[row, 1].legend(loc='upper right', fontsize=fontsize-2)
    for col in range(ncols):
        for m, yr_cut in enumerate(yr_cuts):
            axes[m, col].set_title(yr_cut, fontsize=fontsize)
        axes[row, col].set_ylabel('Pixel Counts', fontsize=fontsize)
        axes[row, col].tick_params(axis='x', labelsize=fontsize-2)
        axes[row, col].tick_params(axis='y', labelsize=fontsize-2)
# set all the ylims
ymin_final, ymax_final = [10**9, 10**9, 10**9], [0, 0, 0]
for i, col in enumerate(range(ncols)):
    for row in range(nrows):
        ymin, ymax = axes[row, col].get_ylim()
        ymin_final[i] = min(ymin_final[i], ymin)
        ymax_final[i] = max(ymax_final[i], ymax)
for i, col in enumerate(range(ncols)):
    for row in range(nrows):
        axes[row, col].set_ylim(ymin_final[i], ymax_final[i])
# set axis labels
axes[m, 0].set_xlabel('Galactic Latitude (deg)', fontsize=fontsize)
axes[m, 1].set_xlabel(r'log$_{10}$ E(B-V)', fontsize=fontsize)
# finalize things
fig.set_size_inches(20, nrows*5)
if save_stuff:
    filename = 'histograms_galLat_ebv_%s_nside%s_%s.png'%(dbname, nside, dither)
    plt.savefig('%s/%s'%(outDir, filename), format= 'png', bbox_inches='tight')
    print('## Saved %s in %s.'%(filename, outDir))
if not dont_show_plots:
    plt.show()
else:
    plt.close('all')

########################################################################################################################
########################################################################################################################
########################################################################################################################
# Finalized cuts
########################################################################################################################
print('#################################################################################################################')
print('#################################################################################################################')
print('#################################################################################################################')
print('#################################################################################################################')
ebv_limit = 0.2
ebv_label = ''
if ebv_cut:
    ebv_label = ' ; EBV<%s'%ebv_limit

# dictionary for the label for the final cut
final_label = {}
for yr in yr_cuts:
    final_label[yr] = 'i>%s%s'%(chosen_cuts[yr], ebv_label)
print('Chosen cuts: %s'%final_label)

# dictionary for the final pixels
final_pixels = {}
for yr in yr_cuts:
    mag_cut = chosen_cuts[yr]
    
    if ebv_cut:
        good_ebv = np.where(ebv_map[iCutPixels[mag_cut][yr]] < ebv_limit)[0]
        final_pixels[yr] = iCutPixels[mag_cut][yr][good_ebv]
        print('%s: throwing away %s pixels further because of EBV cut'%(yr, len(iCutPixels[mag_cut][yr])-len(good_ebv)))
    else:
        final_pixels[yr] = iCutPixels[mag_cut][yr]
        
# print final stats
print('\n#### %s stats: %s: final cuts: %s'%(dbname, dither, final_label))
stats_dict, md_print = calc_stats(bundle=data_bundle, index=final_pixels,
                                  allBandInds=True, return_stuff=True)
md_print = '#### %s stats: %s: final cuts: %s\n%s'%(dbname, dither, final_label, md_print)
################################################################################################
# histogram latitude, extinction
colors = ['m', 'g', 'b', 'r', 'c', 'y']

plt.clf()
nrows, ncols = len(yr_cuts),2
fig, axes = plt.subplots(nrows, ncols)
fig.subplots_adjust(wspace=0.2, hspace=0.3)

max_counts = 0
for i, yr in enumerate(yr_cuts):
    mag_cut = chosen_cuts[yr]
    data_label = 'i>%s%s'%(mag_cut, ebv_label)
    linestyle = 'solid'
    # plot the galactic latitude histogram
    # first plot for no-cut
    lon, lat = hp.pix2ang(nside=nside, ipix=allBandPixels[yr], lonlat=True)
    c = SkyCoord(ra= lon*u.degree, dec=lat*u.degree)
    axes[i, 0].hist(c.galactic.b.degree, label='all-band; all depths>0',
                    bins=bins_b, histtype='step', lw=2, color='k')
    
    # now with the cut
    lon, lat= hp.pix2ang(nside=nside, ipix=final_pixels[yr], lonlat= True)
    c = SkyCoord(ra= lon*u.degree, dec=lat*u.degree)
    axes[i, 0].hist(c.galactic.b.degree, label=data_label,
                    bins=bins_b, histtype='step', lw=2, color=colors[i%len(colors)])

    # ----------
    # print out area for which EBV>0.2, 0.3
    tot = areaPerPixel*len(allBandPixels[yr])
    print('\n%s: no cut\nTotal area (allBand footprint; all depths>0): %.2f deg2'%(yr, tot))
    for lim in [0.2, 0.3]:
        area_here = len(np.where(ebv_map[allBandPixels[yr]]>lim)[0])*areaPerPixel
        print('EBV>%s: Area: %.2f deg2 (%.2f%% of total)'%(lim, area_here, 100.*area_here/tot))
        
    tot = areaPerPixel*len(final_pixels[yr])
    
    print('\n%s:\nTotal area: %s: allBand footprint; all depths>0: %.2f deg2'%(yr, data_label, tot))
    for lim in [0.2, 0.3]:
        area_here = len(np.where(ebv_map[final_pixels[yr]]>lim)[0])*areaPerPixel
        print('EBV>%s: Area: %.2f deg2 (%.2f%% of total)'%(lim, area_here, 100.*area_here/tot))
    # ----------
    # plot the EBV histogram
    # no cut
    cts, _, _ = axes[i, 1].hist(np.log10(ebv_map[allBandPixels[yr]]), label='all-band; all depths>0',
                                bins=bins_ebv, histtype='step', lw=2, color='k')
    max_counts = max(max_counts, max(cts))

    # final cuts
    cts, _, _ = axes[i, 1].hist(np.log10(ebv_map[final_pixels[yr]]),
                                label=data_label,
                                bins=bins_ebv, histtype='step', lw=2,
                                color=colors[i%len(colors)], linestyle=linestyle)
    max_counts = max(max_counts, max(cts))
    
for row in range(nrows):
    x = np.arange(0, max_counts,10)
    for ebv in [0.2, 0.3]: # add lines for constant EBV
        axes[row, 1].plot(np.log10([ebv]*len(x)), x, '-.', label='EBV: %s'%ebv)
    axes[row, 0].legend(loc='upper right', fontsize=fontsize-2)
    axes[row, 1].legend(loc='upper left', fontsize=fontsize-2)
    for col in range(ncols):
        for m, yr_cut in enumerate(yr_cuts):
            axes[m, col].set_title(yr_cut, fontsize=fontsize)
        axes[row, col].set_ylabel('Pixel Counts', fontsize=fontsize)    
        axes[row, col].tick_params(axis='x', labelsize=fontsize-2)
        axes[row, col].tick_params(axis='y', labelsize=fontsize-2)
# set all the ylims
ymin_final, ymax_final = [10**9, 10**9, 10**9], [0, 0, 0]
for i, col in enumerate(range(ncols)):
    for row in range(nrows):
        ymin, ymax = axes[row, col].get_ylim()
        ymin_final[i] = min(ymin_final[i], ymin)
        ymax_final[i] = max(ymax_final[i], ymax)
for i, col in enumerate(range(ncols)):
    for row in range(nrows):
        axes[row, col].set_ylim(ymin_final[i], ymax_final[i])
# set axis labels
axes[m, 0].set_xlabel('Galactic Latitude (deg)', fontsize=fontsize)
axes[m, 1].set_xlabel(r'log$_{10}$ E(B-V)', fontsize=fontsize)
# finalize things
fig.set_size_inches(20, nrows*5)
if save_stuff:
    filename = 'final_footprint_histograms_galLat_ebv_%s_nside%s_%s.png'%(dbname, nside, dither)
    plt.savefig('%s/%s'%(outDir, filename), format= 'png', bbox_inches='tight')
    print('## Saved %s in %s.'%(filename, outDir))
if not dont_show_plots:
    plt.show()
else:
    plt.close('all')

########################################################################################################################
# plot skymaps for each band before and after the depth cut
nTicks = 5
for band in orderBand:
    for yr in yr_cuts:
        mag_cut = chosen_cuts[yr]
        data_label = 'i>%s%s'%(mag_cut, ebv_label)

        plt.clf()
        fig, axes= plt.subplots(1,2)
        temp = copy.deepcopy(data_bundle['%s_%s'%(yr, band)])  # need to copy since will change the mask
        
        # retain the all-band footprint only
        temp.metricValues.mask = True
        temp.metricValues.mask[allBandPixels[yr]] = False
        
        # figure out the color range
        inSurveyIndex = np.where(temp.metricValues.mask == False)[0]
        median = np.median(temp.metricValues.data[inSurveyIndex])
        stddev = np.std(temp.metricValues.data[inSurveyIndex])
        colorMin = median-2.5*stddev
        colorMax = median+2.5*stddev
        increment = (colorMax-colorMin)/float(nTicks)
        ticks = np.arange(colorMin+increment, colorMax, increment)
        
        # plot the no-cut skymap
        plt.axes(axes[0])
        hp.mollview(temp.metricValues.filled(temp.slicer.badval), 
                    flip='astro', rot=(0,0,0) , hold= True,
                    min=colorMin, max=colorMax,
                    title= '%s: No Cut'%yr, cbar=False)
        hp.graticule(dpar=20, dmer=20, verbose=False)
        
        # plot the mag-cut skymap
        temp.metricValues.mask = True
        temp.metricValues.mask[final_pixels[yr]] = False
        
        plt.axes(axes[1])
        hp.mollview(temp.metricValues.filled(temp.slicer.badval), 
                    flip='astro', rot=(0,0,0) , hold= True,
                    min=colorMin, max=colorMax,
                    title= '%s: %s'%(yr, data_label), cbar=False)
        hp.graticule(dpar=20, dmer=20, verbose=False)
            
        # add a color bar
        im = plt.gca().get_images()[0]
        cbaxes = fig.add_axes([0.25, 0.38, 0.5, 0.01]) # [left, bottom, width, height]
        cb = plt.colorbar(im,  orientation='horizontal',
                          ticks=ticks, format='%.2f', cax=cbaxes) 
        cb.set_label('%s-band depth\n(all-band footprint with all depth > 0)'%(band), fontsize=14)
        cb.ax.tick_params(labelsize=14)
        
        fig.set_size_inches(18,18)
        if save_stuff:
            filename = 'final_footprint_skymap_%s_nside%s_%s_%s_%sband.png'%(dbname, nside, dither, yr, band)
            plt.savefig('%s/%s'%(outDir, filename), format= 'png', bbox_inches='tight')
            print('## Saved %s in %s.'%(filename, outDir))
        if not dont_show_plots:
            plt.show()
        else:
            plt.close('all')

if save_stuff:
    # save the final footprint data: area, median/std depth in all bands; 
    # first as a txt file (md table) and then a pickle (dictionary).
    filename = 'footprintStats_%s_%s.txt'%(dbname, dither)
    txt_file = open('%s/%s'%(outDir, filename), 'a')
    txt_file.write(md_print)
    txt_file.close()
    print('## Saved %s in %s.'%(filename, outDir))
    # save pickle too
    filename = 'footprintStats_%s_%s.pickle'%(dbname, dither)
    with open('%s/%s'%(outDir, filename), 'wb') as handle:
        pickle.dump(stats_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    print('## Saved %s in %s.'%(filename, outDir))

    # save/add the final data for the area and the i-band depth for this cadence.
    for yr in yr_cuts:
        header = '| db | cut | Area (deg$^2$) | 5$\sigma$ $i$-band Depth: Median | 5$\sigma$ $i$-band Depth: Std |'
        header += '\n|:--:|:---:|:--------------:|:-------------------------------:|:------------------------------:|'
        db_entry = '\n| %s | %s | %.2f | %.2f | %.2f |'%(dbname, final_label[yr],
                                                         stats_dict['Area (deg$^2$)'][yr],
                                                         stats_dict['5$\sigma$ Depth: Median']['%s_i'%yr],
                                                         stats_dict['5$\sigma$ Depth: Std']['%s_i'%yr]
                                                        )
        # write to the markdown
        filename = 'footprintStats_allcadences_%s_%s.md'%(yr, dither)
        if not os.path.exists('%s/%s'%(outDir_md, filename)):
            to_write = header + db_entry
            update_only = False
        else:
            to_write = db_entry
            update_only = True

        md_file = open('%s/%s'%(outDir_md, filename), 'a')
        md_file.write(to_write)
        md_file.close()

        if update_only:
            print('## Updated %s in %s.'%(filename, outDir_md))
        else:
            print('## Saved %s in %s.'%(filename, outDir_md))

print('\nTime taken: %.2f min'%((time.time()-startTime)/60.))
print('\nAll done.\n')
        
        
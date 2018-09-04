# Goal here is to implement an extinction cut to get rid of the area with high extinction in WFD.
# EBV>0.2 pixels are discarded and, if specified, the final footprint is saved.
#
# Need the coadded depth data for 1yr, 10yr; calculated by running the bas_calc_coadd.sh script.
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
########################################################################################################################
from optparse import OptionParser
parser = OptionParser()
parser.add_option('--nside', dest='nside', type='int',
			help='HEALPix resolution parameter.', default=256)
parser.add_option('--db_path', dest='db_path',
                  help='Path to OpSim db to consider.',
                  default='/global/cscratch1/sd/awan/dbs_wp_unzipped/pontus_2002.db')
parser.add_option('--coadd_data_dir', dest='coadd_data_dir',
                  help='Path to the folder where coadd depth maps are; coadd data with dust extinction and no border masking must exist already.',
                  default='/global/homes/a/awan/LSST/output/coadd_output_noDith/')
parser.add_option('--save_pixels_fIDs',
                  action='store_true', dest='save_pixels_fIDs', default=False,
                  help= 'Set to True if want to save the pixel numbers and the corresponding fieldIDs for the final footprint.')
parser.add_option('--outDir', dest='outDir',
                  help='Path to the folder where all the output should be stored',
                  default='/global/homes/a/awan/LSST/output/')

(options, args) = parser.parse_args()
print('\nOptions: %s'%options)

nside = options.nside
dbpath = options.db_path
data_dir = options.coadd_data_dir
save_pixels_fIDs = options.save_pixels_fIDs
outDir = options.outDir

########################################################################################################################
# set some things up
orderBand = ['u', 'g', 'r', 'i', 'z', 'y']
dbname = dbpath.split('/')[-1].split('.db')[0]

ebv_limit = 0.2

print('\ndata_dir: %s'%data_dir)
print('dbname: %s'%dbname)
print('nside: %s\n'%nside)
print('ebv_limit: %s\n'%ebv_limit)

########################################################################################################################
# get the coadd files and put the data in a bundle
print('## Reading in the data ... \n')
data_bundle = OrderedDict()
for yr_cut in ['1yr', '10yr']:
    if yr_cut=='1yr': file_yearTag = '1yearCut'
    else: file_yearTag = 'fullSurveyPeriod'
        
    for band in orderBand:
        folder = 'coaddM5Analysis_nside%s_withDustExtinction_0pixelRadiusForMasking_%sBand_%s_%s_directory/'%(nside, band,
                                                                                                              dbname,
                                                                                                              file_yearTag)
        path = '%s/%s/unmaskedCoaddData/'%(data_dir, folder)
        filenames = [f for f in os.listdir(path) if f.endswith('.npz')]
        print('Reading %s from\n%s/maskedCoaddData.\n'%(filenames, folder))
        
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

inSurveyIndex = {}
for key in data_bundle.keys():
    inSurveyIndex[key] = np.where(data_bundle[key].metricValues.mask == False)[0]

########################################################################################################################
########################################################################################################################
# calculate some stats and print out a markdown table
areaPerPixel= hp.pixelfunc.nside2pixarea(nside=nside, degrees=True)

def calc_stats(bundle, index, allBandInds=False, return_stuff=False):
    # index must have the same keys as bundle
    if (bundle.keys()!=index.keys()) and not allBandInds:
        raise ValueError('index must have the same keys as bundle:\n%s\n%s'%(bundle.keys(), index.keys()))
        
    if return_stuff and allBandInds: 
        stuff_to_return = {}
        for key in ['5$\sigma$ Depth: Median', '5$\sigma$ Depth: Std', 'Area (deg2)']:
            stuff_to_return[key] = {}
        
    header, sep = '| - ', '| ---- ' 
    med_depth, std_depth, area = '| 5$\sigma$ Depth: Median ', '| 5$\sigma$ Depth: Std ', '| Area (deg2) '
    yr = None
    for key in bundle:
        if yr is None: yr = key.split('yr')[0]+'yr'
        
        current_yr = key.split('yr')[0]+'yr'
        if current_yr!=yr:
            print('%s\n%s\n%s\n%s\n%s\n'%(header, sep, med_depth, std_depth, area))
            header, sep = '| - ', '| ---- '
            med_depth, std_depth, area = '| 5$\sigma$ Depth: Median ', '| 5$\sigma$ Depth: Std ', '| Area (deg2) '
            yr = current_yr
        
        if allBandInds: index_key = current_yr
        else: index_key = key
            
        med = np.nanmedian(bundle[key].metricValues.data[index[index_key]])
        std = np.nanstd(bundle[key].metricValues.data[index[index_key]])
        sarea = (len(index[index_key])*areaPerPixel)
        
        if return_stuff and allBandInds:
            stuff_to_return['5$\sigma$ Depth: Median'][key] = med
            stuff_to_return['5$\sigma$ Depth: Std'][key] = std
            stuff_to_return['Area (deg2)'][index_key] = sarea
            
        header += '| %s '%key
        sep += '| ---- '
        med_depth += '| %.2f '%med
        std_depth += '| %.2f '%std
        area += '| %.2f '%sarea
        
    print('%s\n%s\n%s\n%s\n%s\n'%(header, sep, med_depth, std_depth, area))
    
    if return_stuff: return stuff_to_return
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
# Calculate the stats in the survey region (unmasked; no constraints on depth, i.e., even have negative depths rn).
dat_keys = ['Area (deg2)', '5$\sigma$ Depth: Median', '5$\sigma$ Depth: Std']

########################################################################################################################
# plot galactic latitude and EBV histograms for different cuts
print('\n## Plotting galactic latitude and EBV histograms')
# import EBV map from MAF
opsdb = db.OpsimDatabase(dbpath)
simdata = opsdb.fetchMetricData(['fieldId', 'fieldRA', 'fieldDec', 'night'],
                                sqlconstraint=None)

if dither=='NoDither':
    slicer = slicers.HealpixSlicer(lonCol='fieldRA',
                                   latCol='fieldDec',
                                   latLonDeg=opsdb.raDecInDeg, nside=nside)
else:
    if dither.__contains__('FieldPerNight'):
        s = mafStackers.RandomDitherFieldPerNightStacker(degrees=opsdb.raDecInDeg, randomSeed=1000)
    elif dither.__contains__('PerNight'):
        s = mafStackers.RandomDitherPerNightStacker(degrees=opsdb.raDecInDeg, randomSeed=1000)
    elif dither.__contains__('FieldPerVisit'):
        s = mafStackers.RandomDitherFieldPerVisitStacker(degrees=opsdb.raDecInDeg, randomSeed=1000)
    else:
        raise ValueError('Unsure of what stacker to consdier for %s dithers.'%dither)
    simdata = s.run(simdata)
    dither_timescale = dither.split('Dither')[-1]
    slicer = slicers.HealpixSlicer(lonCol='randomDither%sRa'%dither_timescale,
                                   latCol='randomDither%sDec'%dither_timescale,
                                   latLonDeg=opsdb.raDecInDeg, nside=nside)
# set up dust map to get the ebv map
dustmap = maps.DustMap(nside=nside)
slicer.setupSlicer(simdata)
result = dustmap.run(slicer.slicePoints)
ebv_map = result['ebv']

########################################################################################################################
########################################################################################################################
########################################################################################################################
# Finalized cuts
########################################################################################################################
print('#################################################################################################################')
print('#################################################################################################################')
print('#################################################################################################################')
print('#################################################################################################################')
ebv_label = 'EBV<%s'%ebv_limit
    
final_pixels = {}
for yr in ['1yr', '10yr']:
    good_ebv = np.where(ebv_map[allBandPixels[yr]] <= ebv_limit)[0]
    final_pixels[yr] = allBandPixels[yr][good_ebv]
    
# print final stats
print('\n#### %s stats: %s: final cuts: %s'%(dbname, dither, ebv_label))
calc_stats(bundle=data_bundle, index=final_pixels, allBandInds=True)

################################################################################################
# histogram latitude, extinction
bins_b = np.arange(-90, 90, 0.5)
bins_ebv = np.arange(-3.5, 0.5, 0.05)
colors = ['m', 'g', 'b', 'r', 'c', 'y']
fontsize = 14
plt.clf()
fig, axes = plt.subplots(2,2)
fig.subplots_adjust(wspace=0.2, hspace=0.3)
max_counts = 0
for i, yr in enumerate(['1yr', '10yr']):
    data_label = '%s'%(ebv_label)

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
    
for row in [0, 1]:
    x = np.arange(0,max_counts,10)
    for ebv in [0.2, 0.3]: # add lines for constant EBV
        axes[row, 1].plot(np.log10([ebv]*len(x)), x, '-.', label='EBV: %s'%ebv)
    axes[row, 0].set_xlabel('Galactic Latitude (deg)', fontsize=fontsize)
    axes[row, 1].set_xlabel(r'log$_{10}$ E(B-V)', fontsize=fontsize) 
    axes[row, 1].set_ylim(0,max_counts)
    
    for col in [0, 1]:
        axes[0, col].set_title('1yr', fontsize=fontsize)
        axes[1, col].set_title('10yr', fontsize=fontsize)
        axes[row, col].legend(loc='upper right', fontsize=fontsize-2)
        axes[row, col].set_ylabel('Pixel Counts', fontsize=fontsize)    
        axes[row, col].tick_params(axis='x', labelsize=fontsize-2)
        axes[row, col].tick_params(axis='y', labelsize=fontsize-2)
plt.suptitle('allBand footprint; all depths>0', fontsize=fontsize)
fig.set_size_inches(20,10)
plt.show()

################################################################################################
# depth histogram
bins_depth_full = np.arange(0, 27, 0.1)
bins_depth_cut = bins_depth_full[bins_depth_full>20]
plt.clf()
fig, axes = plt.subplots(2,2)
fig.subplots_adjust(wspace=0.2, hspace=0.3)

max_counts = 0
for i, yr in enumerate(['1yr', '10yr']):
    data_label = '%s'%(ebv_label)

    linestyle = 'solid'
    for j, band in enumerate(orderBand):
        # first plot for no-cut
        axes[i, 0].hist(data_bundle['%s_%s'%(yr, band)].metricValues.data[allBandPixels[yr]],
                     label='%s-band'%band, bins=bins_depth_full,
                     histtype='step', lw=2, color=colors[j])
        # now with cut
        axes[i, 1].hist(data_bundle['%s_%s'%(yr, band)].metricValues.data[final_pixels[yr]],
                     label='%s-band'%band, bins=bins_depth_cut,
                     histtype='step', lw=2, color=colors[j])
for row in [0, 1]:
    axes[row, 0].set_xlabel('5$\sigma$ Coadded depth', fontsize=fontsize)
    for col in [0, 1]:
        axes[row, col].legend(loc='upper left', fontsize=fontsize-2)
        axes[row, col].set_ylabel('Pixel Counts', fontsize=fontsize)    
        axes[row, col].tick_params(axis='x', labelsize=fontsize-2)
        axes[row, col].tick_params(axis='y', labelsize=fontsize-2)

axes[0, 0].set_title('1yr: all-depths>0', fontsize=fontsize)
axes[0, 1].set_title('1yr: %s'%data_label, fontsize=fontsize)
axes[1, 0].set_title('10yr: all-depths>0', fontsize=fontsize)
axes[1, 1].set_title('10yr: %s'%data_label, fontsize=fontsize)

plt.suptitle('allBand footprint; all depths>0', fontsize=fontsize)
fig.set_size_inches(20,10)
plt.show()

########################################################################################################################
# plot skymaps for each band before and after the depth cut
nTicks = 5
for band in orderBand:
    for yr in ['1yr', '10yr']:
        data_label = '%s'%(ebv_label)
        
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
        colorMin = median-1.5*stddev
        colorMax = median+1.5*stddev
        increment = (colorMax-colorMin)/float(nTicks)
        ticks = np.arange(colorMin+increment, colorMax, increment)
        
        # plot the no-cut skymap
        plt.axes(axes[0])
        hp.mollview(temp.metricValues.filled(temp.slicer.badval), 
                    flip='astro', rot=(0,0,0) , hold= True,
                    min=colorMin, max=colorMax,
                    title= '%s: all-band; all depths>0'%yr, cbar=False)
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
        plt.show()
        
if save_pixels_fIDs:
    import pickle
    for yr in ['10yr']:
        data_to_save = {}
        data_to_save['pixNum'] = final_pixels[yr]
    
        # now find all the fieldIDs that have observations
        fID_list = []
        for pixel in data_to_save['pixNum']:
            indObsInPixel = slicer._sliceSimData(pixel)
            fID_list += list(simdata[indObsInPixel['idxs']]['fieldId']) # fieldIDs corresponding to pixel
        
        data_to_save['fieldIDs'] = np.unique(fID_list)
        filename = '%sFootprint_%s_nside%s_%s_ebv<%s.pickle'%(yr, dbname, nside, dither, ebv_limit)
            
        with open('%s/%s'%(outDir, filename), 'wb') as handle:
            pickle.dump(data_to_save, handle, protocol=pickle.HIGHEST_PROTOCOL)
            
        print('\n## Saved %s in %s.\n'%(filename, outDir))
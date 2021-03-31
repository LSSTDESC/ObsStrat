##############################################################################################################
# This script calculate the overlap area between eg-footprint and specified survey footprint.
# (based on the notebook used to calculate overlaps for 2018 white paper call.)
##############################################################################################################
import numpy as np
import os
import healpy as hp
import pandas as pd
import time
import matplotlib.pyplot as plt
#######################################################################################
from optparse import OptionParser
parser = OptionParser()
parser.add_option('--outdir', dest='outdir',
                  help='Path to the folder where all the output should be stored; \
                  if the folder doesnt exist, it will be created.')
parser.add_option('--eg-path', dest='eg_path',
                  help='Path to the folder with all the (LSST) eg-footprints.')
parser.add_option('--nside-lsst', dest='nside_lsst', type='int',
                  help='HEALPix resolution parameter for the LSST footprint.')
parser.add_option('--other-footprint-path', dest='footprint_path',
                  help='Path to the file containing the footprint with which overlap is to be calculated.')
parser.add_option('--survey', dest='survey_name',
                  help='Survey name for which the footprint is inputted; the format \
                  is assumed to be a csv file with one column containing the \
                  pixel numbers that fall in the survey footprint.')
parser.add_option('--nside-other', dest='nside_other', type='int',
                  help='HEALPix resolution parameter for the non-LSST footprint.')
parser.add_option('--dbs-order-path', dest='dbs_order_path',
                  help='Path to the file with the order of LSST sims dbs to consider.')
#######################################################################################
start_time = time.time()
# get the inputs
(options, args) = parser.parse_args()
print(options)
outdir = options.outdir
eg_path = options.eg_path
nside_lsst = options.nside_lsst
footprint_path = options.footprint_path
survey_name = options.survey_name
nside_other = options.nside_other
dbs_order_path = options.dbs_order_path

# set up directory for output data
os.makedirs(outdir, exist_ok=True)

# read in the other-footprint pixels
other_footprint_pixels = pd.read_csv(footprint_path).values.flatten()
if len(other_footprint_pixels) <= 0:
    raise ValueError('something is wrong with the footprint path for %s' % survey_name)
# --------------------------------------------------------------------------------
# helper functions
# --------------------------------------------------------------------------------
def get_area(nside, npixels):
    return hp.nside2pixarea(nside=nside, degrees=True) * npixels

def plot_skymap(tag, hpixel_in, nside, return_masked_arr=False):
    # set up a masked array
    npix = hp.nside2npix(nside=nside)
    vals = np.zeros(npix) + 5000
    footprint = vals.view(np.ma.MaskedArray)
    footprint.mask = [True] * npix
    footprint.fill_value = np.nan
    footprint.mask[hpixel_in] = False

    area = get_area(nside=nside, npixels=len(hpixel_in) )
    # plot things out
    plt.clf()
    hp.mollview(footprint, flip='astro', rot=(0,0,0), min=0, max=10,
                title='%s footprint: %.2f deg2; nside %s' % (tag, area, nside)
               )
    hp.graticule(dpar=20, dmer=20, verbose=False)
    plt.savefig('%s/skymap_%s_nside%s.png' % (outdir, tag, nside))
    plt.close('all')

    if return_masked_arr:
        return footprint
# --------------------------------------------------------------------------------
# lets plot the skymap out
masked_arr = plot_skymap(tag=survey_name, hpixel_in=other_footprint_pixels,
                         nside=nside_other, return_masked_arr=True)

# see if need to change the map resolution
if nside_other != nside_lsst:
    masked_arr.fill_value = -100 # to enable udgrad
    # now change the resolution
    masked_arr = hp.ud_grade(masked_arr, nside_out=nside_lsst)
    other_footprint_pixels = np.where(masked_arr.data > 0)[0] # i.e. ignore pixels with 
    # lets plot out the skymap again
    plot_skymap(tag=survey_name,
                hpixel_in=other_footprint_pixels, nside=nside_lsst)
    print('## converted %s footprint from nside %s to nside %s' % (survey_name, nside_other, nside_lsst))
masked_arr = []

# set up all the years to consider
yr_cuts = [1, 3, 6, 10]
# get the dbs order
dbs_in_order = pd.read_csv(dbs_order_path).values.flatten()
ndbs = len(dbs_in_order)
print('## %s dbs to consider.' % ndbs)
# calculate the area per pixel
area_per_pix = hp.nside2pixarea(nside=nside_lsst, degrees=True)

# the order in which we want to save things is: db-name, Y1 overlap, Y3, Y6, Y10
overlaps = {}
for db in dbs_in_order:
    print('## running things for %s ...' % db)
    if db.__contains__('.db'): db = db.split('.db')[0]
    for yr in yr_cuts:
        # locate the file
        fname = [f for f in os.listdir(eg_path) if f.endswith('txt') \
                 and f.__contains__(db) and f.__contains__('yr%s_' % yr)]
        if len(fname) != 1:
            raise ValueError('somethings wrong - expecting 1 file but got %s' % len(fname))
        fname = fname[0]
        # load the lsst eg-footprint
        eg_mask_path = '%s/%s' % (eg_path, fname)
        eg_mask = np.array( np.genfromtxt(eg_mask_path) )
        # find the indices for the in-survey pixels
        lsst_pix = np.where(eg_mask == 0)[0]

        # now calculate the overlap
        overlap = list(set(lsst_pix)-(set(lsst_pix)-set(other_footprint_pixels)))
        # add area to the right list
        if yr not in overlaps:
            # initiate
            overlaps[yr] = []
        # add the latest overlap
        overlaps[yr].append( float( '%.2f' % (area_per_pix * len(overlap))) )

# set up the structure to save everything
df_dict = {'dbname': dbs_in_order}
for yr in overlaps:
    df_dict['Y%s overlap with %s (deg2)' % (yr, survey_name)] = overlaps[yr]
df = pd.DataFrame(df_dict)

# save the data
yr_tag = ''
for yr in yr_cuts: yr_tag += 'Y%s_' % yr
fname = 'overlaps-%s_%s%sdbs.csv' % (survey_name, yr_tag, ndbs)
df.to_csv('%s/%s' % (outdir, fname), index=False)
print('## saved %s' % fname)

# all done.
print('## time taken for %s overlaps: %.2f min\n' % (survey_name, (time.time() - start_time )/60 ) )
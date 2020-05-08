##############################################################################################################
# This script runs compare_versions to compare baseline between v1.3 and v1.4.
#
# Humna Awan: humna.awan@rutgers.edu
#
##############################################################################################################
import time
import os

from helper_compare_versions_misc_metrics import compare_versions

# set up
start_time = time.time()
nside = 64 # not the most high resolution but its faster..

yr_cut = 10
outdir = '/global/homes/a/awan/LSST/lsstRepos/ObsStrat/postwp/results-plots+/compare-big-sky-dust_v1.4-v1.5/'
os.makedirs(outdir, exist_ok=True)

dbpath_dict = {}
dbpath_dict['v1.4'] = '/global/cscratch1/sd/awan/dbs_post_wp_v3/footprint_big_sky_dustv1.4_10yrs.db'
dbpath_dict['v1.5'] = '/global/cscratch1/sd/awan/dbs_post_wp_orig_v4/footprints/footprint_big_sky_dustv1.5_10yrs.db'

dbname = 'big-sky-dust'

reference_version = 'v1.4'
order_of_versions = ['v1.4', 'v1.5']

ilims = {'v1.4': 25.9,  'v1.5': 25.9}

ebvlims = 0.2

# run the code
compare_versions(outdir=outdir, dbpath_dict=dbpath_dict, dbname=dbname,
                 reference_version=reference_version, order_of_versions=order_of_versions,
                 nside=nside, yr_cut=yr_cut, ilims=ilims, ebvlims=ebvlims)

print('## time taken: %.f min' % ((time.time() - start_time) / 60.))
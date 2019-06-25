##############################################################################################################
# This scripts saves the rotTelPos, rotSkyPos, PA for the new FBS output with rot dithers.
#
# Humna Awan: humna.awan@rutgers.edu
#
##############################################################################################################
import os
import pandas as pd
import lsst.sims.maf.db as db
import lsst.sims.maf.stackers as stackers
import time

import matplotlib
matplotlib.use('Agg')
import os
import pandas as pd
import time
import matplotlib.pyplot as plt
import pylab
import numpy as np

fontsize = 18
pylab.rcParams['axes.labelsize'] = fontsize
pylab.rcParams['xtick.labelsize'] = fontsize-2
pylab.rcParams['ytick.labelsize'] = fontsize-2
pylab.rcParams['legend.fontsize'] = fontsize
pylab.rcParams['axes.linewidth'] = 2
#######################################################################################
from optparse import OptionParser
parser = OptionParser()
parser.add_option('--dbfile', dest='dbfile',
                  help='Path to the folder with the db to consider.')
parser.add_option('--outdir', dest='outdir',
                  help='Path to the folder with the db to consider.',
                  default='/global/homes/a/awan/awan/desc/rot_output/')

(options, args) = parser.parse_args()
dbfile = options.dbfile
outdir = options.outdir
#######################################################################################
rot_rand_seed = 42

dbname = dbfile.split('/')[-1].split('.db')[0]
print('Running save_rot_data for %s'%(dbname))
time0 = time.time()
# ----------------------------------------------------------
# connect to the database
opsdb = db.OpsimDatabase(dbfile)
# ----------------------------------------------------------
# WFD only
# fetch the data: columns need for the rotational dither and parallactic angle stacker
colnames=['fieldRA', 'fieldDec', 'fieldId', 'rotTelPos', 'rotSkyPos', \
         'observationStartMJD', 'observationStartLST', 'note']
simdata = opsdb.fetchMetricData(colnames=colnames, sqlconstraint=None)
# add PA
s = stackers.ParallacticAngleStacker(degrees=opsdb.raDecInDeg)
simdata = s.run(simdata)
# add fieldID; adds a column named opsimFieldId
s = stackers.OpSimFieldStacker(degrees=opsdb.raDecInDeg)
simdata = s.run(simdata)

# get the WFD visits
nobs = len(simdata['note'])
wfd_ind = [f for f in range(nobs) if not simdata['note'][f].__contains__('DD')]

# wrap rotTelPos
ind = simdata['rotTelPos'] > 100
simdata['rotTelPos'][ind] -= 360

# ----------------------------------------------------------
data_dict = {}
for col in ['fieldRA', 'fieldDec', 'fieldId', 'rotTelPos', \
            'rotSkyPos', 'PA', 'note']:
    if col == 'fieldId':
        data_dict[col] = simdata['opsimFieldId'][wfd_ind]
    else:
        data_dict[col] = simdata[col][wfd_ind]

print('min, max rotTelPos: %s, %s'%(min(data_dict['rotTelPos']), max(data_dict['rotTelPos']) ) )
# ----------------------------------------------------------
# save the output array
filename = '%s_data.csv'%(dbname)
df = pd.DataFrame(data_dict)
df.to_csv('%s/%s'%(outdir, filename), index=False)
print('\nSaved %s'%filename)

print('Total time taken: %.2f min'%((time.time()-time0)/60.))

#######################################################################################
time0 = time.time()

# now read in the data
db = filename.split(filename)[0]
simdata = data_dict
# ----------------------------------------------------------
# set up tp plot the histograms for the data
# set up bins
bins_telpos = np.arange(-95, 95, 1)
bins_skypos = np.arange(0, 360, 1)

# plot the histogram for the data
nrows, ncols = 2, 1
fig, axes = plt.subplots(nrows=nrows, ncols=ncols)
plt.subplots_adjust(hspace=0.3, wspace=0.1)

fids = np.unique(simdata['fieldId'])
nfids = len(fids)
colors = [(0, np.random.random(), 0) for i in range(nfids)]
# loop over all the fields
for j, fid in enumerate(fids):
    # --------------------------------------
    ind = np.where(simdata['fieldId'] == fid)[0]
    axes[0].hist(simdata['rotTelPos'][ind], bins=bins_telpos, color=colors[j],
                    histtype='step', lw=2)
    # undithered rotSkyPos
    axes[1].hist(simdata['rotSkyPos'][ind],
                    bins=bins_skypos, color=colors[j],
                    histtype='step', lw=2)
# plot details
axes[0,].set_xlabel('rotTelPos (degrees)')
axes[1].set_xlabel('rotSkyPos (degrees)')
for row in range(nrows):
    axes[row].set_ylim([0, 50])
    axes[row].grid()
    axes[row].set_ylabel('Counts')
plt.suptitle(r'$\bf{%s}: %s \ fields$'%(dbname.replace('_', ' '), len(fids)), fontsize=fontsize)
fig.set_size_inches(10, 10)
filename = 'compare_rot_angles_%s_perfield.png'%(dbname)
plt.savefig('%s/%s'%(outdir, filename), format= 'png', bbox_inches='tight')
print('Saved %s'%filename)
plt.show()
print('Time taken for the plot: %.2f min\n'%((time.time()-time0)/60.))
##############################################################################################################
# This scripts plots the distributions in different fields for dithered and undithered rotTelPos and
# rotSkyPos for different dbs. The code was initially run in an iPython notebook but it takes a really long
# time to plots things.
#
# Humna Awan: humna.awan@rutgers.edu
#
##############################################################################################################
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
parser.add_option('--dbname', dest='dbname',
                  help='db name.',
                  default='baseline2018a')
parser.add_option('--outdir', dest='outdir',
                  help='Path to the folder with the db to consider.',
                  default='/global/homes/a/awan/desc/wp_descDithers_csvs/compare_rot_dith/')

(options, args) = parser.parse_args()
dbname = options.dbname
outdir = options.outdir
#######################################################################################
time0 = time.time()
# get all the files
data_files = [f for f in os.listdir(outdir) if f.endswith('csv') and f.__contains__(dbname)]
if len(data_files)>1:
    raise ValueError('Found more than one csv file: %s'%data_files)

file = data_files[0]
# now read in the data
print(file)
db = file.split('_data.csv')[0]
simdata = pd.read_csv('%s/%s'%(outdir, file))
print('\n## Read in data for %s'%(dbname))

# ----------------------------------------------------------
# set up tp plot the histograms for the data
dith_skypos = np.array(simdata['randomDitherPerFilterChangeRotTelPos']-simdata['PA'])
# wrap the dithered rotSkyPos
ind = dith_skypos < 0
dith_skypos[ind] += 360
simdata['randomDitherPerFilterChangeRotSkyPos'] = dith_skypos.copy()

# set up bins
bins_telpos = np.arange(-90, 90, 1)
bins_skypos = np.arange(0, 360, 1)
print(dbname)

# plot the histogram for the data
nrows, ncols = 2,2
fig, axes = plt.subplots(nrows=nrows, ncols=ncols)
plt.subplots_adjust(hspace=0.3, wspace=0.1)

fids = np.unique(simdata['fieldId'])
nfids = len(fids)
colors = [(0, np.random.random(), 0) for i in range(nfids)]
# loop over all the fields
for j, fid in enumerate(fids):
    # --------------------------------------
    ind = np.where(simdata['fieldId'] == fid)[0]
    axes[0, 0].hist(simdata['rotTelPos'][ind], bins=bins_telpos, color=colors[j],
                    histtype='step', lw=2)
    # dithered rotTelPos
    axes[0, 1].hist(simdata['randomDitherPerFilterChangeRotTelPos'][ind],
                    #label=db,
                    color=colors[j],
                    bins=bins_telpos, histtype='step', lw=2)
    # undithered rotSkyPos
    axes[1, 0].hist(simdata['rotSkyPos'][ind],
                    bins=bins_skypos, color=colors[j],
                    histtype='step', lw=2)
    # dithered rotSkyPos
    axes[1, 1].hist(simdata['randomDitherPerFilterChangeRotSkyPos'][ind],
                    bins=bins_skypos, color=colors[j],
                    histtype='step', lw=2)
# plot details
ymax_tel, ymax_sky = 220, 50
for row in range(nrows):
    for col in range(ncols):
        axes[0, col].set_ylim(0, ymax_tel)
        axes[1, col].set_ylim(0, ymax_sky)
        #axes[row, col].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        axes[row, 0].set_ylabel('Counts')
        axes[0, col].set_xlabel('rotTelPos (degrees)')
        axes[1, col].set_xlabel('rotSkyPos (degrees)')
        axes[row, col].grid()
#axes[0, 1].legend(bbox_to_anchor=(1.,1))
axes[0, 0].set_title(r'$\bf{undithered}$', fontsize=fontsize, y=1.08)
axes[0, 1].set_title(r'$\bf{perFilter \ rotational \ dithers}$', fontsize=fontsize, y=1.08)
plt.suptitle(r'$\bf{%s}: %s \ fields$'%(dbname.replace('_', ' '), len(fids)), fontsize=fontsize)
fig.set_size_inches(20,10)
filename = 'compare_rot_angles_%s_perfield_nodith_wdith.png'%(dbname)
plt.savefig('%s/%s'%(outdir, filename), format= 'png', bbox_inches='tight')
print('Saved %s'%filename)
plt.show()
print('Time taken for the plot: %.2f min\n'%((time.time()-time0)/60.))
##############################################################################################################
# this script reads in the overlap data and creates a plot.
##############################################################################################################
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import time
from matplotlib.lines import Line2D
# -----------------------------------------------------------------------------------------------------
import sys
sys.path.append('/global/homes/a/awan/my-utils/')
from pyutils import settings
# -----------------------------------------------------------------------------------------------------
#######################################################################################
from optparse import OptionParser
parser = OptionParser()
parser.add_option('--fbs-version', dest='fbs_version',
                  help='Path to the folder the overlap csv files are stored.')
#######################################################################################
start_time = time.time()
# get the inputs
(options, args) = parser.parse_args()
print(options)
fbs_version = options.fbs_version

if fbs_version == 'v1.5':
    from helper_dbcategories import db_catmap_v15 as db_catmap
elif fbs_version == 'v1.6':
    from helper_dbcategories import db_catmap_v16 as db_catmap
elif fbs_version == 'v1.7':
    from helper_dbcategories import db_catmap_v17 as db_catmap
else:
    raise ValueError('dont have folder map for the specified fbs-version: %s' % fbs_version)
    
datadir = '/global/homes/a/awan/LSST/lsstRepos/ObsStrat/postwp/paper-data/overlaps_%s/' % fbs_version
# set up the surveys
surveys = ['4MOST', 'DESI' ,'euclid']
# yrs to consider
yrs = [1, 3, 6, 10]

# read in the data
datas = {}
for survey in surveys:
    fname = [f for f in os.listdir(datadir) if f.endswith('csv') and f.__contains__('%s_' % survey)][0]
    print('## reading in %s' % fname)
    datas[survey] = pd.read_csv('%s/%s' % (datadir, fname))
    ndbs = int(fname.split('dbs')[0].split('_')[-1])
print('## %s dbs to consider' % ndbs)

# set up the xlabels
xlabels = datas[survey]['dbname'].values
for i, xlabel in enumerate(xlabels):
    if xlabel.endswith('.db'): xlabels[i] = xlabel.split('.db')[0]

# set up the colors for the db catagories
colors_dbs = ['indianred', 'mediumslateblue', 'olive', 'orangered',
             'turquoise', 'brown', 'goldenrod', 'dodgerblue', 'darkorchid',
             'y', 'palevioletred', 'teal', 'sandybrown', 'forestgreen',
              '#9467bd', '#8c564b', '#e377c2', '#17becf'
            ]
grp_colors = []
custom_lines, grp_labels = [], []
for db in xlabels:
    db_found = False
    for i, group in enumerate(db_catmap.keys()):
        if db in db_catmap[group]:
            db_found = True
            grp_colors += [colors_dbs[i]]
            if group not in grp_labels:
                grp_labels += [group]
                custom_lines += [Line2D([0], [0], color=colors_dbs[i], lw=10) ]
    if db_found is False:
        grp_colors += ['#7f7f7f']
        grp = 'no-grp'
        if grp not in grp_labels:
            grp_labels += [grp]
            custom_lines += [Line2D([0], [0], color='#7f7f7f', lw=10) ]

if ndbs != len(grp_colors):
    raise ValueError('something is wrong - missing a group!')
    
colors = ['#2ca02c', '#ff7f0e', '#1f77b4', '#d62728',]
#colors = ['m', 'b', 'g', 'k']
shapes = ['d', 'o', 's', 'v']

# plot and group things by folder-group
plt.clf()
nrows, ncols = 1, len(surveys)
fig, axes = plt.subplots(nrows, ncols)
plt.subplots_adjust(wspace=0.05, hspace=0.2, top=0.9)

# loop over the survey
for j, survey in enumerate(surveys):
    # loop over the yr
    for i, yr in enumerate(yrs):
        col = 'Y%s overlap with %s (deg2)' % (yr, survey)
        axes[j].plot( datas[survey][col] , range(ndbs), '%s' % shapes[i], color=colors[i], label='Y%s' % yr)
    # set title
    axes[j].set_title('overlap with %s' % survey)
# add legend for the points
axes[-1].legend(bbox_to_anchor=(1.,1))
# plot details
for ncol in range(ncols):
    axes[ncol].set_yticks( range(ndbs) )
    axes[ncol].set_xlabel( r'deg$^2$' )
axes[0].set_yticklabels(xlabels )
for col in range(1, ncols):
    axes[col].set_yticklabels([])

# now add legend for the db catagories
# color-code the yticks
for ytick, color in zip(axes[0].get_yticklabels(), grp_colors):
    ytick.set_color(color) 
legend = axes[0].legend(custom_lines, grp_labels, bbox_to_anchor=(1.8, 1.08),
                        frameon=True, ncol=6, title='db groups')
plt.setp(legend.get_title(), fontsize=16)
# fig size
plt.gcf().set_size_inches(20, 20 * (ndbs/50) )
# save figure
filename = '%s_compare-overlaps_%sdbs_grouped.png' % (fbs_version, ndbs)
plt.savefig('%s/%s' % (datadir, filename), format='png', bbox_inches='tight')
print('\n# Saved %s' % filename)
plt.close('all')

print('## Time taken: %.2f min\n' % ((time.time() - start_time )/60 ) )
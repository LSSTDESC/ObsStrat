##############################################################################################################
# This script runs the egFootprintMetric with the option to get the i-band depth map, calculates eg footprint
# area and depth statistics, and saves them. Also saves the depth data bundle.
#
# Humna Awan: humna.awan@rutgers.edu
#
##############################################################################################################
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import os
import pandas as pd
from helper import folder_map
from matplotlib.lines import Line2D
# -----------------------------------------------------------------------------------------------------
fontsize = 18
rcparams = {}
rcparams['figure.figsize'] = (10, 6)
rcparams['axes.labelsize'] = fontsize
rcparams['legend.fontsize'] = fontsize-4
rcparams['axes.titlesize'] = fontsize
rcparams['axes.linewidth'] = 2
rcparams['axes.grid'] = True
for axis in ['x', 'y']:
    rcparams['%stick.labelsize'%axis] = fontsize-2
    rcparams['%stick.direction'%axis] = 'in'
    rcparams['%stick.major.size'%axis] = 5.5
    rcparams['%stick.minor.size'%axis] =  3.5
    rcparams['%stick.major.width'%axis] = 2
    rcparams['%stick.minor.width'%axis] = 1.5
rcparams['xtick.top'] = True
rcparams['ytick.right'] = True
for key in rcparams: mpl.rcParams[key] = rcparams[key]
# -----------------------------------------------------------------------------------------------------
outdir = '/global/homes/a/awan/LSST/lsstRepos/ObsStrat/code/postwp/'
data_dir = '/global/cscratch1/sd/awan/lsst_output/post_wp_output/summary_data/'

# set up for plots
colors = ['m', 'b', 'g', 'k']
shapes = ['o', 'v', 's', 'd']

colors_cm = ['indianred', 'mediumslateblue', 'olive', 'orangered', 'black',
             'turquoise', 'brown', 'goldenrod', 'dodgerblue', 'darkorchid', 'y', 'palevioletred', 'teal', 'sandybrown']
# -----------------------------------------------------------------------------------------------------
# set up directory for bundle data
results_dir = '%s/results/' % outdir
os.makedirs(results_dir, exist_ok=True)
# -----------------------------------------------------------------------------------------------------
# read in the data
data = {}
yr_label = {}
for yr in [1, 3, 6, 10]:
    files = [f for f in os.listdir( data_dir ) if f.endswith('csv') and f.__contains__('y%s_' % yr) and f.startswith('eg_') ]
    for file in files:
        print( 'Readig in %s' % file )
        key = 'yr%s_%s' % (yr, file.split('_')[4])
        data[ key ] = pd.read_csv('%s/%s' % (data_dir, file))
        yr_label[ key ] =  r'Y%s (i$>$%s) ' % (yr, file.split('_')[4].split('limi')[-1])
# add ngal data if available for all dbs
for yr in [1, 3, 6, 10]:
    files = [f for f in os.listdir( data_dir ) if f.endswith('csv') and f.__contains__('y%s_' % yr) and f.startswith('ngal_') ]
    for file in files:
        print( 'Readig in %s' % file )
        key = 'yr%s_%s' % (yr, file.split('_')[6])
        temp = pd.read_csv('%s/%s' % (data_dir, file))
        if len(temp['dbname'].values) == len(data[key]['dbname'].values):            
            data[key] = pd.merge(temp, data[key], left_on='dbname', right_on='dbname',how='outer')

# make sure the headers make sense
for key in data:
    data[key] = data[key].rename(columns=lambda x: x.strip())
        
# order data in a specific way
from itertools import chain
order_group = list (chain(*folder_map.values()) )
df_b = pd.DataFrame({'dbname' : order_group})

for key in data:
    data[key] = pd.merge(df_b, data[key], left_on='dbname',right_on='dbname',how='outer')
# set up the colors
folder_colors = []
custom_lines, grp_labels = [], []
for i, group in enumerate(folder_map.keys()):
    for key in folder_map[group]:
        folder_colors += [colors_cm[i]]
    grp_labels += [group]         
    custom_lines += [Line2D([0], [0], color=colors_cm[i], lw=10) ]
# -----------------------------------------------------------------------------------------------------
# plot and group things by folder-group
xlabels = np.array( data[list(data.keys())[0]]['dbname'].values )
ndbs_tot = len(xlabels)
for j, to_plot in enumerate( ['Area (deg2)', '$i$-band depth: median', '$i$-band depth: std', 'Ngal'] ):
    plt.clf()
    nrows, ncols = 1, 2
    fig, axes = plt.subplots(nrows, ncols)
    plt.subplots_adjust(wspace=1., hspace=0.2, top=0.9)
    for i, yr in enumerate( data ):
        if to_plot in data[yr]:
            data_to_plot = data[yr][to_plot].values
            if i + j == 0:
                ind = range( len( data_to_plot ) )
                mid = int(len(ind)/2)
                ndbs = ndbs_tot
            # plot
            axes[0].plot(data_to_plot[ ind[0:mid] ] , range(ndbs)[0:mid], '%s'%shapes[i], color=colors[i], label=yr_label[yr]) 
            axes[1].plot(data_to_plot[ ind[mid:] ] , range(ndbs)[mid:], '%s'%shapes[i], color=colors[i], label=yr_label[yr]) 
    # plot 18K line
    if to_plot == 'Area (deg2)':
        axes[0].plot(( [18000] * ndbs)[0:mid], range(ndbs)[0:mid], 'r--', lw=2, label='18K') 
        axes[1].plot(( [18000] * ndbs)[mid:], range(ndbs)[mid:], 'r--', lw=2, label='18K') 
        
    # plot details
    axes[-1].legend(bbox_to_anchor=(1.,1))
    axes[0].set_yticks( range(ndbs)[0:mid] )
    axes[0].set_yticklabels(xlabels[ ind[0:mid] ] )
    axes[1].set_yticks( range(ndbs)[mid:] )
    axes[1].set_yticklabels(xlabels[ ind[mid:] ] ) #rotation=90)
    # color-code the yticks
    for ytick, color in zip(axes[0].get_yticklabels(), folder_colors[0:mid]):
        ytick.set_color(color)    
    for ytick, color in zip(axes[1].get_yticklabels(), folder_colors[mid:]):
        ytick.set_color(color)
    # set up the ylims
    min_val, max_val = axes[0].get_xlim()
    other_min, other_max = axes[1].get_xlim()
    min_val = min([min_val, other_min])
    max_val = max([max_val, other_max])
    for i in range(2):
        axes[i].set_xlabel(r'%s' % to_plot)
        axes[i].set_xlim(min_val, max_val)
    # set up the legend
    axes[0].legend(custom_lines, grp_labels, bbox_to_anchor=(2.8, 1.1), frameon=True, ncol=7)
    # figure size
    plt.gcf().set_size_inches(20, 15)
    # set up to save fig
    if to_plot == 'Area (deg2)':
        plot_label = 'area'
    if to_plot == '$i$-band depth: median':
        plot_label = 'iband_med-depth'
    if to_plot == '$i$-band depth: std':
        plot_label = 'iband_std-depth'
    if to_plot == 'Ngal':
        plot_label = 'ngal'
    # save fig
    filename = 'compare_%s_%sdbs_grouped.png'%(plot_label, ndbs_tot)
    plt.savefig('%s/%s'%(results_dir, filename), format= 'png', bbox_inches='tight')
    print('\n# Saved %s'%filename)
    plt.close('all')
# -----------------------------------------------------------------------------------------------------
# plot and group things yr1 area
xlabels = np.array( data[list(data.keys())[0]]['dbname'].values )
ndbs = len(xlabels)
for j, to_plot in enumerate( ['Area (deg2)', '$i$-band depth: median', '$i$-band depth: std', 'Ngal'] ):
    plt.clf()
    nrows, ncols = 1, 2
    fig, axes = plt.subplots(nrows, ncols)
    plt.subplots_adjust(wspace=1.0, hspace=0.2, top=0.9)
    
    for i, yr in enumerate( data ):
        if to_plot in data[yr]:
            data_to_plot = data[yr][to_plot].values
            if i + j == 0:
                print('Ordering by ascending %s for %s' % (to_plot, yr) )
                ind = np.argsort( data_to_plot )
                mid = int(len(ind)/2)
            # plot  
            axes[0].plot(data_to_plot[ ind[0:mid] ] , range(ndbs)[0:mid], '%s'%shapes[i], color=colors[i], label=yr_label[yr]) 
            axes[1].plot(data_to_plot[ ind[mid:] ] , range(ndbs)[mid:], '%s'%shapes[i], color=colors[i], label=yr_label[yr]) 
    # plot 18K limit for area
    if to_plot == 'Area (deg2)':
        axes[0].plot(( [18000] * ndbs)[0:mid], range(ndbs)[0:mid], 'r--', lw=2, label='18K') 
        axes[1].plot(( [18000] * ndbs)[mid:], range(ndbs)[mid:], 'r--', lw=2, label='18K') 
    # plot details
    axes[-1].legend(bbox_to_anchor=(1.,1))
    axes[0].set_yticks( range(ndbs)[0:mid] )
    axes[0].set_yticklabels(xlabels[ ind[0:mid] ] )
    axes[1].set_yticks( range(ndbs)[mid:] )
    axes[1].set_yticklabels(xlabels[ ind[mid:] ] ) #rotation=90)
    # axis lims
    min_val, max_val = axes[0].get_xlim()
    other_min, other_max = axes[1].get_xlim()
    min_val = min([min_val, other_min])
    max_val = max([max_val, other_max])
    for i in range(2):
        axes[i].set_xlabel(r'%s' % to_plot)
        axes[i].set_xlim(min_val, max_val)
    # color the ticks by folder name
    folder_colors = np.array(folder_colors)
    for ytick, color in zip(axes[0].get_yticklabels(), folder_colors[ ind[0:mid] ]):
        ytick.set_color(color)
    for ytick, color in zip(axes[1].get_yticklabels(), folder_colors[ ind[mid:] ]):
        ytick.set_color(color)
    # overall legend
    axes[0].legend(custom_lines, grp_labels, bbox_to_anchor=(2.8, 1.1), frameon=True, ncol=7)
    # fig size
    plt.gcf().set_size_inches(20, 15)
    # set up to save
    if to_plot == 'Area (deg2)':
        plot_label = 'area'
    if to_plot == '$i$-band depth: median':
        plot_label = 'iband_med-depth'
    if to_plot == '$i$-band depth: std':
        plot_label = 'iband_std-depth'
    if to_plot == 'Ngal':
        plot_label = 'ngal'
    # save fig
    filename = 'compare_%s_%sdbs_ordered_by_area.png'%(plot_label, ndbs_tot)
    plt.savefig('%s/%s'%(results_dir, filename), format= 'png', bbox_inches='tight')
    print('\n# Saved %s'%filename)
    plt.close('all')

# -----------------------------------------------------------------------------------------------------
outdir = '%s/stats/' % results_dir
os.makedirs(outdir, exist_ok=True)
# save the data
cuts = {}
for key in data:
    yr_cut = key.split('_')[0].split('yr')[-1]
    mag_cut = key.split('_')[1].split('limi')[-1]
    cuts[ '%syr' % yr_cut ] = 'i>%s ; EBV<0.2' % mag_cut

    areas = data[key]['Area (deg2)']
    med =  data[key]['$i$-band depth: median']
    std =  data[key]['$i$-band depth: std']
    # store what to write
    to_write = ''
    for i, db in enumerate( data[key]['dbname'].values ):
        to_write += '| %s | %s | %.2f | %.2f |\n' % (db, cuts[ '%syr' % yr_cut ], areas[i], med[i])
    # now write to the file
    fname = 'stats_Y%s.txt' % (yr_cut)
    txt_file = open('%s/%s' % (outdir, fname), 'a')
    txt_file.write(to_write)
    txt_file.close()

    print('Saved %s' % fname)
# Goal here is to automate saving out dithered positions (RA, Dec, rotTelPos) for non-MAF users.
#
# We calculate the translational and rotational dithers for various cadences and
# save the output as a csv file.  These dithers are largely the same as in DC1/DC2:
#   - Translational dithers:
#       - WFD: large random offsets (as large as 1.75 deg) applied after every visit.
#       - DD: small random offsets (as large as 7 arcmin) applied after every visit.
#       - Else: no dithers, so `fieldRA`, `fieldDec` are returned.
#   - Rotational dithers:
#       - All surveys (WFD, DD, else): random between -90, 90 degrees applied after every filter change.
#                                      (Break from DC2: Some visits dont get dithered since they are
#                                                        forced outside the rotator range.
#                                                        See RotStacker info for details.)
# Supports OpSim V3/V4 outputs.
#
# Saved file format:
#   .csv file with four columns: obsIDcol, 'descDitheredRA', 'descDitheredDec', 'descDitheredRotTelPos'
#    where
#        obsIDcol = 'observationId' for V4 outputs and 'obsHistID' for V3 outputs.
#
#    Saved filename = descDithers_<database name>.csv
#
# Humna Awan: humna.awan@rutgers.edu
#
########################################################################################################################
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import lsst.sims.maf
import lsst.sims.maf.db as db
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.metricBundles as metricBundles
import lsst.sims.maf.stackers as stackers
import time
import pandas as pd
import os
import numpy as np
import datetime
########################################################################################################################
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-q', '--quiet',
                  action= 'store_false', dest= 'print_progress', default=True,
                  help= 'Do not print stuff out')
parser.add_option('--dbs_path', dest='dbs_path',
                  help='Path to the directory that contains the .db files; could have non-.db files.',
                  default='/global/cscratch1/sd/awan/dbs_wp_unzipped')
parser.add_option('--outDir', dest='outDir',
                  help='Path to the folder where all the output should be stored; directory should exist already.',
                  default='/global/homes/a/awan/desc/wp_descDithers_csvs')
parser.add_option('--db_files_only', dest='db_files_only',
                  help='List of names of the db files to run. E.g. <baseline2018a.db, pontuns_2002.db>',
                  default=None)
parser.add_option('--rot_rand_seed', dest='rot_rand_seed', type='int',
                  help='Seed for random number generator for rotational dithers.', default=42)
parser.add_option('--trans_rand_seed', dest='trans_rand_seed', type='int',
                  help='Seed for random number generator for translational dithers.', default=42)
parser.add_option('--save_plots',
                  action='store_true', dest='save_plots', default=False,
                  help= 'Use to save the histogram for descDithers in outDir.')
parser.add_option('--show_diagnostic_plots',
                  action='store_true', dest='show_diagnostic_plots', default=False,
                  help= 'Use to show histogram of added dithers.')
parser.add_option('--compress_csvs',
                  action='store_true', dest='compress_csvs', default=False,
                  help= 'Use to compress the csv file outputs.')

########################################################################################################################
startTime = time.time()
(options, args) = parser.parse_args()
print('\nOptions: %s'%options)

# read in the inputs
print_progress = options.print_progress
dbs_path = options.dbs_path
outDir = options.outDir
db_files_only = options.db_files_only
rot_rand_seed = options.rot_rand_seed
trans_rand_seed = options.trans_rand_seed
save_plots = options.save_plots
show_diagnostic_plots = options.show_diagnostic_plots
compress_csvs = options.compress_csvs

if db_files_only is not None:
    db_files_only = [f.strip() for f in list(db_files_only.split(','))]
########################################################################################################################
startTime_0 = time.time()
readme = '#########################################################################################\n'
readme = '%s\n'%(datetime.date.isoformat(datetime.date.today()))
readme += 'Running with lsst.sims.maf.__version__: %s\n\n'%lsst.sims.maf.__version__
readme += '## save_csv_dithers run:\n\ndbs_path= %s\n'%dbs_path
readme += 'outDir: %s\n'%outDir
readme += 'db_files_only: %s\n'%db_files_only
readme += 'rot_rand_seed: %s\ntrans_rand_seed: %s\n'%(rot_rand_seed, trans_rand_seed)
readme += 'print_progress: %s\nshow_diagnostic_plots: %s\n'%(print_progress, show_diagnostic_plots)
readme += 'compress_csvs: %s\n'%compress_csvs

dbfiles = [f for f in os.listdir(dbs_path) if f.endswith('db')]  # select db files
if print_progress: print('Found files: %s\n'%dbfiles)

if db_files_only is not None:
    dbfiles = [f for f in dbfiles if f in db_files_only]  # select db files

readme += '\nReading for files: %s\n\n'%dbfiles
if print_progress and db_files_only is not None: print('Running over: %s\n'%dbfiles)

for i, dbfile in enumerate(dbfiles): # loop over all the db files
    startTime = time.time()
    if (i!=0): readme = ''
    readme += '%s'%dbfile

    if print_progress: print('Starting: %s\n'%dbfile)

    opsdb = db.OpsimDatabase('%s/%s'%(dbs_path, dbfile)) # connect to the database

    # specify the column names to get from the db file
    colnames = ['proposalId', 'observationId', 'fieldRA', 'fieldDec', 'rotTelPos']
    propIDcol, obsIDcol= 'proposalId', 'observationId'

    if (opsdb.opsimVersion=='V3'):
        # V3 outputs have somewhat different column names
        colnames = ['propID', 'obsHistID', 'fieldRA', 'fieldDec', 'rotTelPos']
        propIDcol, obsIDcol= 'propID', 'obsHistID'

    # get the data
    simdata = opsdb.fetchMetricData(colnames=colnames, sqlconstraint=None)

    # set up to run the stackers that add columns for translational and rotational dithers.
    metric = metrics.PassMetric()  # want to access the database; no analysis needed
    slicer = slicers.OneDSlicer(sliceColName='night', binsize=1, verbose=print_progress)   # essentially accessing all nights
    sqlconstraint = None

    resultsDb = db.ResultsDb(outDir=outDir)
    ########################################################################################################################
    # set up metric bundle to run stackers for large translational dithers + rotational dithers
    if print_progress: print('Setting up for WFD translational dithers + rot dithers.')
    bgroup= {}
    stackerList = [stackers.RandomDitherFieldPerVisitStacker(degrees=opsdb.raDecInDeg,
                                                            randomSeed=trans_rand_seed),
                   stackers.RandomRotDitherPerFilterChangeStacker(degrees=opsdb.raDecInDeg,
                                                                  randomSeed=rot_rand_seed)]

    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint=sqlconstraint,
                                        stackerList=stackerList)
    bgroup['WFD'] = metricBundles.MetricBundleGroup({0: bundle}, opsdb, outDir=outDir,
                                                    resultsDb=resultsDb, saveEarly=False,
                                                    verbose=print_progress)
    # run the bundle
    bgroup['WFD'].runAll()

    # set up the bundle for small translational dithers
    if print_progress: print('\nSetting up for DD translational dithers.')
    chipSize= 1.75*2/15
    chipMaxDither= chipSize/2.
    stackerList = [stackers.RandomDitherFieldPerVisitStacker(maxDither= chipMaxDither,
                                                             degrees=opsdb.raDecInDeg,
                                                             randomSeed=trans_rand_seed)]
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint=sqlconstraint,
                                         stackerList=stackerList)

    bgroup['DD'] = metricBundles.MetricBundleGroup({0: bundle}, opsdb, outDir=outDir,
                                                   resultsDb=resultsDb, saveEarly=False,
                                                   verbose=print_progress)
    # run the bundle
    bgroup['DD'].runAll()

    ########################################################################################################################
    # access the relevant columns
    dithered_RA, dithered_Dec = {}, {}
    for key in bgroup:
        dithered_RA[key] = bgroup[key].simData['randomDitherFieldPerVisitRa']
        dithered_Dec[key] = bgroup[key].simData['randomDitherFieldPerVisitDec']

    dithered_rotTelPos = bgroup['WFD'].simData['randomDitherPerFilterChangeRotTelPos']

    ########################################################################################################################
    # diagnostic plots
    if show_diagnostic_plots:
        # histograms of dithers
        fig, axes = plt.subplots(nrows=1, ncols=3)

        for key in bgroup:
            # ra
            axes[0].hist(dithered_RA[key]-simdata['fieldRA'],
                         label='%s dithers: delRA'%key, histtype='step', lw=2, bins= 30)

            # dec
            axes[1].hist(dithered_Dec[key]-simdata['fieldDec'],
                         label='%s dithers: delDec'%key, histtype='step', lw=2)

        # tel pos
        axes[2].hist(dithered_rotTelPos-simdata['rotTelPos'],
                     label='rot dithers: rotTelPos', histtype='step', lw=2)
        for ax in axes:
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            ax.set_ylabel('Counts')

        axes[0].legend()
        axes[1].legend()

        if opsdb.raDecInDeg: unitlabel = 'degrees'
        else: unitlabel = 'radians'

        axes[0].set_xlabel('delRA (%s)'%unitlabel)
        axes[1].set_xlabel('delDec (%s)'%unitlabel)
        axes[2].set_xlabel('delRotTelPos (%s)'%unitlabel)

        plt.title(dbfile)
        fig.set_size_inches(20,5)

    ########################################################################################################################
    # initiate the final arrays as undithered fieldRA, fieldDec as nonWFD, nonDDF should remain unchanged
    descDitheredRA = simdata['fieldRA'].copy()
    descDitheredDec = simdata['fieldDec'].copy()
    descDitheredRot = simdata['rotTelPos'].copy()

    # need to find the indices for WFD vs. DD observations since we are adding different
    # translational dithers for WFD/DDF visits + none for other surveys
    propIds, propTags = opsdb.fetchPropInfo()
    # ok work with WFD visits now
    ind_WFD = np.where(simdata[propIDcol]==propTags['WFD'])[0]
    if print_progress:
        tot= len(simdata)
        print('Total visits: ', tot)
        print('propTags: ', propTags)
        print('%s WFD visits out of total %s'%(len(ind_WFD), tot))

    descDitheredRA[ind_WFD] = dithered_RA['WFD'][ind_WFD]
    descDitheredDec[ind_WFD] = dithered_Dec['WFD'][ind_WFD]

    # work with DD visits now
    ind_DD = np.where(simdata[propIDcol]==propTags['DD'])[0]
    if print_progress:
        print('%s DD visits out of total %s'%(len(ind_DD), tot))

    descDitheredRA[ind_DD] = dithered_RA['DD'][ind_DD]
    descDitheredDec[ind_DD] = dithered_Dec['DD'][ind_DD]

    # add rotational dithers to everything
    descDitheredRot = dithered_rotTelPos

    ########################################################################################################################
    # diagnostic plots
    if show_diagnostic_plots or save_plots:
        # histograms of desc dithered positions
        fig, axes = plt.subplots(nrows=1, ncols=3)

        _, bins, _ = axes[0].hist(descDitheredRA, label='descDitheredRA', histtype='step', lw=2)
        axes[0].hist(simdata['fieldRA'], label='fieldRA', histtype='step', lw=2, bins=bins)

        _, bins, _ = axes[1].hist(descDitheredDec, label='descDitheredDec', histtype='step', lw=2)
        axes[1].hist(simdata['fieldDec'], label='fieldDec', histtype='step', lw=2, bins=bins)

        _, bins, _ = axes[2].hist(descDitheredRot, label='descDitheredRot', histtype='step', lw=2)
        axes[2].hist(simdata['rotTelPos'], label='rotTelPos', histtype='step', lw=2, bins=bins)

        if opsdb.raDecInDeg: xlabel = 'degrees'
        else: xlabel = 'radians'

        for ax in axes:
            ax.legend()
            ax.set_xlabel(xlabel)
            ax.set_ylabel('Counts')

        plt.suptitle(dbfile)
        fig.set_size_inches(20,5)

        if save_plots:
            filename = 'hist_descDithers_%s.png'%(dbfile.split('.db')[0])
            plt.savefig('%s/%s'%(outDir, filename), format='png', bbox_inches='tight')
            readme += '\nSaved hist for descDithers in %s.'%filename

            if print_progress:
                print('\nSaved hist plot in %s'%filename)

        if show_diagnostic_plots:
            plt.show()
        else:
            plt.close('all')

    ########################################################################################################################
    # save the columns as a csv file.
    d= {obsIDcol: simdata[obsIDcol],
        'descDitheredRA': descDitheredRA, 'descDitheredDec': descDitheredDec,
        'descDitheredRotTelPos': descDitheredRot}

    filename = 'descDithers_%s.csv'%(dbfile.split('.db')[0])
    if compress_csvs:
        filename = '%s.gz'%(filename)
        pd.DataFrame(d).to_csv('%s/%s'%(outDir, filename), index=False, compression='gzip')
    else:
        pd.DataFrame(d).to_csv('%s/%s'%(outDir, filename), index=False)

    readme += '\nSaved the dithers in %s'%filename
    readme += '\nTime taken: %.2f (min)\n\n'%((time.time()-startTime)/60.)

    if print_progress:
        print('\nSaved the dithers in %s'%filename)
        print('Time taken: %.2f (min)\n\n'%((time.time()-startTime)/60.))

    readme_file = open('%s/readme.txt'%(outDir), 'a')
    readme_file.write(readme)
    readme_file.close()

# mark the end in the readme.
readme_file= open('%s/readme.txt'%(outDir), 'a')
readme_file.write('All done. Total time taken: %.2f (min)\n\n'%((time.time()-startTime_0)/60.))
readme_file.close()
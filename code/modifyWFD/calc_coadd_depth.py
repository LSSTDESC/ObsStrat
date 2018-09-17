# Goal here is to calculate the 5sigma coadded depth for all the cadences (or just the baseline and pontuns2002).
# Uses the coaddM5Analysis code from mafContrib. Run for all bands. No border masking.
#
# Humna Awan: humna.awan@rutgers.edu
#
########################################################################################################################
import matplotlib
matplotlib.use('Agg')
import os
import time
from mafContrib import coaddM5Analysis
########################################################################################################################
from optparse import OptionParser
parser = OptionParser()
parser.add_option('--nside', dest='nside', type='int',
			help='HEALPix resolution parameter.', default=256)
parser.add_option('--yr_cut', dest='yr_cut', type='int',
			help='Year cut.', default=None)
parser.add_option('--no_dither',
                  action='store_true', dest='no_dither', default=False,
                  help= 'Run the analyis without any translational dithers.')
parser.add_option('--noMWdust',
                  action='store_true', dest='noMWdust', default=False,
                  help= 'Use the tag to not include MW dust extinction.')
parser.add_option('--baseline_and_wide_only',
                  action='store_true', dest='baseline_and_wide_only', default=False,
                  help='Use the tag to consider only baseline2018a and pontus2002.')
parser.add_option('--dbs_path', dest='dbs_path',
                  help='Path to the folder with the unzipped OpSim dbs.',
                  default='/global/cscratch1/sd/awan/dbs_wp_unzipped')
parser.add_option('--outDir', dest='outDir',
                  help='Path to the folder where all the output should be stored. The folder must already exist.',
                  default='/global/homes/a/awan/LSST/output/coadd_output_allwps/')
parser.add_option('--specific_db', dest='specific_db',
                  help='Specific db to run the analysis for.',
                  default=None)
parser.add_option('--slair',
                  action='store_true', dest='slair', default=False,
                  help= 'Use the tag to specify a slair run.')

(options, args) = parser.parse_args()
nside = options.nside
yr_cut = options.yr_cut
no_dither = options.no_dither
dbs_path = options.dbs_path
outDir = options.outDir
includeDustExtinction = not options.noMWdust
baseline_and_wide_only = options.baseline_and_wide_only
specific_db = options.specific_db
slair = options.slair

########################################################################################################################
# set some things up
bands = ['u', 'g', 'r', 'i', 'z', 'y']

if yr_cut==10: yr_cut=None

dbfiles = [f for f in os.listdir(dbs_path) if f.endswith('db')]
if baseline_and_wide_only:
    dbfiles = [f for f in dbfiles if ((f=='baseline2018a.db') or (f=='pontus_2002.db'))]
    print('Running over %s'%dbfiles)

if specific_db is not None:
    if specific_db in dbfiles:
        dbfiles = [specific_db]
    else:
        raise ValueError('%s is not in the dbs_path directory.'%specific_db)

if no_dither: dith = 'NoDither'
else: dith = 'RandomDitherPerNight'
    
startTime_0 = time.time()
# no border masking.
for dbfile in dbfiles:
    runName = dbfile.split('.db')[0]
    print('runName: %s'%runName)
    for filterBand in bands:
        startTime = time.time()
        print('filter= %s\n'%filterBand)
        out = coaddM5Analysis.coaddM5Analysis(path=outDir, dbfile='%s/%s'%(dbs_path, dbfile),
                                              runName=runName, WFDandDDFs=False, slair=slair,
                                              specifiedDith = dith,
                                              nside=nside, filterBand=filterBand,
                                              includeDustExtinction=includeDustExtinction, saveunMaskedCoaddData=True,
                                              pixelRadiusForMasking=0, cutOffYear=yr_cut,
                                              plotSkymap=True, plotCartview=False,
                                              unmaskedColorMin=None, unmaskedColorMax=None,
                                              maskedColorMin=None, maskedColorMax=None,
                                              nTicks=5, plotPowerSpectrum=True,
                                              showPlots=False, saveFigs=True,
                                              almAnalysis=False)
        print('#### Total taken for %s: %s-band: %.2f min'%(runName, filterBand, (time.time()-startTime)/60.))
        
print('\nTotal time: %.2f min'%((time.time()-startTime_0)/60.))
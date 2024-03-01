# Collecting Routines for use in MAF analysis
# Code written by:
# R. Mandelbaum
# R. Hlozek
# Rubin Sim Notebooks contributors https://github.com/lsst/rubin_sim_notebooks/graphs/contributors
# A. Kannawadi

import os
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import pandas as pd
import sqlite3
import rubin_sim
import rubin_sim.maf as maf
from rubin_sim.scheduler.utils import SkyAreaGenerator
from rubin_sim.data import get_baseline
import scipy.special as sc

def load_maf_map(fname,nside=64):
    surveyAreas = SkyAreaGenerator(nside=nside)
    map_footprints, map_labels = surveyAreas.return_maps()
    fin = maf.MetricBundle(
        metric=maf.ExgalM5(),
        slicer=maf.HealpixSubsetSlicer(nside=nside, use_cache=False, hpid=np.where(map_labels == "lowdust")[0]),
        constraint='note not like "DD%" and note not like "twilight_near_sun"'
    )
    fin._setup_metric_values()
    #counts/map values
    opmap=fin.metric_values
    mafmap=np.copy(opmap.data)
    mafmap[opmap.mask]=0
    mafmap[mafmap>0]=1
    
    return mafmap

def maf_maps_to_fits(fname_in, fname_out,nside=64):
    usemap = load_maf_map(fname_in,nside=nside)
    # save this into fits file:
    hp.write_map(fname_out, usemap, overwrite=True)
    print("Written: %s"%fname_out)


# Here we define a function for some of the metric plots we want to show.
# TO DO: need to rename this function since it doesn't actually always plot the metrics
def metric_plots(use_run_name, use_opsim_fname, use_metric=maf.ExgalM5(), use_color_min=None, use_color_max=None,
                year=10,nside=64, use_filter="i", return_map=False):
    # use_run_name says which OpSim DB we want to use, e.g. `baseline_v2.1_10yrs` - will also be used for labels
    # use_opsim_fname says where it lives, e.g. `/global/cfs/cdirs/lsst/groups/CO/rubin_sim/sim_baseline/baseline_v2.1_10yrs.db`
    surveyAreas = SkyAreaGenerator(nside=nside)
    map_footprints, map_labels = surveyAreas.return_maps()
    days = year*365.25
    # Here the constraint on use of i-band data, exclusion of DDFs, time limitations, and avoiding twilight exposures 
    constraint_str='filter="YY" and note not like "DD%" and night <= XX and note not like "twilight_near_sun" '
    constraint_str = constraint_str.replace('XX','%d'%days)
    constraint_str = constraint_str.replace('YY','%s'%use_filter)
    
    # Just some optional plotting stuff
    if use_color_min is not None and use_color_max is not None:
        plot_dict={"color_min": use_color_min, "color_max": use_color_max, "x_min": use_color_min, "x_max": use_color_max}
    else:
        plot_dict=None
        
    # Decide what summary statistics we want
    my_summary_stats = [maf.MedianMetric(), maf.MeanMetric(), maf.RmsMetric(), maf.PercentileMetric(percentile=25), maf.PercentileMetric(percentile=75)]
    
    # First, define a MetricBundle object.
    depth_map_bundle = maf.MetricBundle(
        metric=use_metric,
        # Exclude the galactic plane
        slicer=maf.HealpixSubsetSlicer(nside=nside, use_cache=False, hpid=np.where(map_labels == "lowdust")[0]),
        constraint=constraint_str,
        run_name=use_run_name,
        summary_metrics=my_summary_stats,
        plot_dict=plot_dict
    )
    
    bd = maf.metricBundles.make_bundles_dict_from_list([depth_map_bundle])
    bgroup = maf.MetricBundleGroup(
        bd, use_opsim_fname
    )
    bgroup.run_all()
    
    if return_map:
        map= depth_map_bundle.metric_values
        output_map = np.copy(map.data)
        output_map[map.mask] = 0
        return bgroup, bd, output_map
    else:
        return bgroup, bd

def coeff_solve(ilim,meanz):

    from scipy.optimize import curve_fit
    ((a, b), _) = curve_fit(lambda x,a,b:3*a*x+3*b, ilim, meanz)
    return a, b

def p(z, ilim, a=0.044444444444444446, b=-0.7644444444444444):
    z0 = a*ilim + b
    return 0.5/z0*(z/z0)*(z/z0)*np.exp(-z/z0)



def percentile_bin_boundaries(*, ilim, num_bins=5, a=0.044444444444444446, b=-0.7644444444444444, **kwargs):
    """Return the tomographic bin boundaries.
    
    The number of returned values will be one more than ``num_bins``.
    The first entry should always be 0 and the last entry should be np.inf.
    """
    z0 = a*ilim + b
    return sc.gammaincinv(3, np.arange(0, 1 + 0.1/num_bins, 1./num_bins))*z0

def fixed_bin_boundaries(*, num_bins=5, zmin=0, zmax=2, a=0.044444444444444446, b=-0.7644444444444444, **kwargs):
    return np.linspace(zmin, zmax, num_bins+1)


def mean_z(ilim, num_bins=5,a=0.044444444444444446, b=-0.7644444444444444, tomographic_bin_boundaries=percentile_bin_boundaries, **kwargs):
    """Compute the mean redshift for each tomographic bin
    given the limiting i-band magnitude.
    
    This function uses the already calibrated value for z_0,
    but can be overridden by passing in the a and b coefficients,
    optionally.
    
    Parameters
    ----------
    ilim: `float`
        The limiting magnitude in i-band
    num_bins: `int`
        Number of tomographic bins
    a: `float`, optional
        The coefficient in z0 = a*ilim + b
    b: `float`, optional
        The offset in z0 = a*ilim + b
    tomographic_bin_boundaries: `Callable`
        A function that gives the edges of the tomographic bins.
        Additional arguments to pass to ``tomographic_bin_boundaries``
        can be passed as keyword arguments.
        
    Returns
    -------
    mean_z: `list` [`float`]
        A list of ``num_bins`` in length with mean redshift in each tomographic bin defined
        by ``tomographic_bin_boundaries``.
    """
    z0 = a*ilim + b
    t_edges = tomographic_bin_boundaries(ilim=ilim, num_bins=num_bins, a=a, b=b, **kwargs)/z0
    t_bins = zip(t_edges[:-1], t_edges[1:])
    return [3*z0*(sc.gammainc(4, t2)-sc.gammainc(4, t1))/(sc.gammainc(3, t2)-sc.gammainc(3, t1)) for (t1,t2) in t_bins]

def sensitivity(num_bins, fiducial_ilim=24.1, a=0.044444444444444446, b=-0.7644444444444444, tomographic_bin_boundaries=percentile_bin_boundaries, **kwargs):
    """Sensitivity of mean redshift to limiting magnitude (d<z>/dilim)
    
    Parameters
    ----------
    num_bins: `int`
        Number of tomographic bins
    fiducial_ilim: `float`, optional
        The fiducial limiting magnitude at which to calculate the sensitivity.
    a: `float`, optional
        The coefficient in z0 = a*ilim + b
    b: `float`, optional
        The offset in z0 = a*ilim + b
    tomographic_bin_boundaries: `Callable`
        A function that gives the edges of the tomographic bins.
        Additional arguments to pass to ``tomographic_bin_boundaries``
        can be passed as keyword arguments.
        
    Returns
    -------
    sensitivity: `list` [`float`]
        The derivative of mean_z wrt ilim for each tomographic bin defined
        by ``tomographic_bin_boundaries``.
    """
    z0 = a*fiducial_ilim + b
    t_edges = tomographic_bin_boundaries(num_bins=num_bins, ilim=fiducial_ilim, **kwargs)/z0
    t_bins = zip(t_edges[:-1], t_edges[1:])
    return [3*a*(sc.gammainc(4, t2)-sc.gammainc(4, t1))/(sc.gammainc(3, t2)-sc.gammainc(3, t1)) for t1,t2 in t_bins]

def my_total_power_metric(map, ell_max=30):
    mean_val = np.mean(map[map>0])
    map[map>0] -= mean_val
        
    cl = hp.anafast(map)
    ell = np.arange(np.size(cl))
    return np.sum((2*ell[ell<ell_max]+1)*cl[ell<ell_max])

# A utility to plot summary stats for strategies as a function of year, given a dataframe from the above routines.
def plot_metric_by_year(df, stat_name,years=None, y_axis_label=None,ylog=False):

    if years!=None:
        print(years)
        year_vals = years
        print(year_vals)
    else:
        year_vals = np.array(list(set(df['Year'])))
        print(year_vals)

    strategies=list(set(df['Strategy']))
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Offset axes:
    offsets = 0.05*np.arange(0,len(strategies))
    offsets -= np.mean(offsets)
    offset_index = 0

    # y-axis label handling
    if y_axis_label is None:
        y_axis_label = stat_name

        
    ## put in line style stuff
    for s in strategies:
        #yvals = np.array([float(val) for val in df[stat_name][df['Strategy']==s]])
        yvals = np.array(df[stat_name][df['Strategy']==s])
        ax.plot(year_vals+offsets[offset_index], yvals, label=s)
        offset_index += 1
    plt.xlabel('Year')
    plt.ylabel(y_axis_label)
    if ylog: plt.yscale('log')
    plt.legend()
    plt.show()
    
# A utility to plot summary stats for strategies as a function of year, given a dataframe from the above routines.
def plot_meanz_metrics_by_year_ak(df, years, num_bins=5,y_axis_label=None):
    year_vals = years
    strategies=list(set(df['Strategy']))
    # fig = plt.figure()
    # ax = fig.add_subplot(111)

    # Offset axes:
    offsets = 0.01*np.arange(0,len(strategies))
    offsets -= np.mean(offsets)
    offset_index = 0

    # y-axis label handling
    if y_axis_label is None:
        y_axis_label = 'Mean z'

    ls = ['--','-',':','-.','-']
    lw = [1,2,2,2,1]
        
    ## put in line style stuff
    cols = ['lightcoral', 'mediumpurple','deepskyblue','teal','forestgreen','r','g','k','c','m']
    offset = 0.05
    fig, axs = plt.subplots(3,len(year_vals),sharex=True, figsize=(15,12))
  
    print(np.shape(axs))
    for sy,year in enumerate(year_vals):
        axs[0][sy].set_title(f'Year {year}', fontsize=15)
        #axs[1][sy].set_xlabel('Tomographic Bin', fontsize=15)
        axs[2][sy].set_xlabel('Tomographic Bin', fontsize=15)
        axs[2][sy].set_ylabel('Cl bias', fontsize=15)

        for scount,s in enumerate(strategies):
            for bin in range(num_bins):      
                meanz=df['Mean z bin'][df['Strategy']==s][df['Year']==year].values[0][bin]
                stdz=df['Std z bin'][df['Strategy']==s][df['Year']==year].values[0][bin]
                if (sy==len(year_vals)-1 and bin==0):
                    axs[0][sy].plot(bin+offset*scount,meanz, marker='*',label=s,color=cols[scount])
                    axs[1][sy].plot(bin+offset*scount,stdz, marker='*',label=s,color=cols[scount])
                else:
                    axs[0][sy].plot(bin+offset*scount,meanz, marker='*',color=cols[scount])
                    axs[1][sy].plot(bin+offset*scount,stdz, marker='*',color=cols[scount])
        
        nbins_use = np.shape(df['Used meanz'][df['Strategy']==s][df['Year']==year].values[0])[0]
        #print(nbins_use, 'nbins_use')
        for binn in range(nbins_use):
            clbias = df['Clbias'][df['Strategy']==s][df['Year']==year].values[0][binn]
        
            if (sy==len(year_vals)-1 and bin==0):
                axs[2][sy].plot(binn+offset*scount,clbias, marker='*',label=s,color=cols[scount])
            else:
                axs[2][sy].plot(binn+offset*scount,clbias, marker='*',color=cols[scount])

    #     
        #axs[sy].set_yscale('log')

    
    axs[0][0].set_ylabel('Mean z', fontsize=15)
    axs[1][0].set_ylabel('Std z', fontsize=15)
    axs[1][sy].legend(loc='upper left',fontsize=10)
    #axs[0].show()
# A utility to plot summary stats for strategies as a function of year, given a dataframe from the above routines.

# First define a routine to run across a list of years and produce a dataframe
def get_year_by_year_metrics_ak(year_list, name_list, sim_list, use_filter="i"):
    overall_names = []
    overall_years = []
    overall_meds = []
    overall_means = []
    overall_std = []
    overall_iqr = []
    overall_meanzbins=[]
    overall_stdzbins=[]
    overall_clbias = []
    meanz_usecl= []
    for year in year_list:
        for i in range(len(sim_list)):
            bgroup, bd = metric_plots(name_list[i], sim_list[i], year=year, use_filter=use_filter)
            overall_names.append(name_list[i]) # strategy name
            overall_years.append(year) 
            overall_meds.append(bd[list(bd.keys())[0]].summary_values['Median']) # median i-band mags
            overall_means.append(bd[list(bd.keys())[0]].summary_values['Mean'])   # mean i-band mags        
            overall_std.append(bd[list(bd.keys())[0]].summary_values['Rms']) # rms of the i-band mags
            overall_iqr.append(bd[list(bd.keys())[0]].summary_values['75th%ile']-bd[list(bd.keys())[0]].summary_values['25th%ile']) 
            meanz = mean_z(ilim=float(bd[list(bd.keys())[0]].summary_values['Mean']), num_bins=5)
            overall_meanzbins.append(meanz) # mean z in each tomographic bin
            stdz = [float(sens)*float(bd[list(bd.keys())[0]].summary_values['Rms']) for sens in sensitivity(fiducial_ilim=float(bd[list(bd.keys())[0]].summary_values['Mean']), num_bins=5)]
            clbias, meanz_use = compute_Clbias(meanz,stdz)
            overall_clbias.append(clbias)
            meanz_usecl.append(meanz_use)
            # We send the i-band magnitude - 1 (so one mag brighter than the output of the i-band limiting magnitude) to Arun's sensitivity code
            # we then multiply the sensitivity from Arun's code by the std of the i-band mag to get the std of the z in each bin
            # print(bd[list(bd.keys())[0]].summary_values['Rms'], 'rms')
            # print([float(sens) for sens in sensitivity(fiducial_ilim=float(bd[list(bd.keys())[0]].summary_values['Mean']), num_bins=5)])
            # print(stdz)
            # testing
            #stdz = [float(sens) for sens in sensitivity(fiducial_ilim=float(bd[list(bd.keys())[0]].summary_values['Mean']), num_bins=5)]
            overall_stdzbins.append(stdz)
            
    df = pd.DataFrame(list(zip(overall_names, overall_years, overall_meds, overall_means, overall_std, overall_iqr, overall_meanzbins,overall_stdzbins, overall_clbias, meanz_usecl)), 
                  columns=['Strategy', 'Year', 'Median i-band depth', 'Mean i-band depth', 'Std i-band depth', 'IQR i-band depth', 'Mean z bin', 'Std z bin','Clbias','Used meanz'])
    return df

def get_year_by_year_metrics(year_list, name_list, sim_list):
    overall_names = []
    overall_years = []
    overall_meds = []
    overall_means = []
    overall_std = []
    overall_iqr = []
    overall_meanzbins=[]
    overall_stdzbins=[]
    overall_clbias = []
    meanz_usecl= []
    filter_list=["u","g","r","i","z","y"]

    for year in year_list:
        for count,i in enumerate(range(len(sim_list))):
            overall_names.append(name_list[i]) # strategy name
            overall_years.append(year)
            for use_filter in filter_list:
                bgroup, bd = metric_plots(name_list[i], sim_list[i], year=year, use_filter=use_filter)
 
                mag = bd[list(bd.keys())[0]].summary_values['Mean']
            overall_meds.append(bd[list(bd.keys())[0]].summary_values['Median']) # median i-band mags
            overall_means.append(imag)   # mean i-band mags        
            overall_std.append(bd[list(bd.keys())[0]].summary_values['Rms']) # rms of the i-band mags
            overall_iqr.append(bd[list(bd.keys())[0]].summary_values['75th%ile']-bd[list(bd.keys())[0]].summary_values['25th%ile']) 
           
            zbins=1
            dzdminterp, meanzinterp=compute_dzfromdm(zbins, imag,year, 'SRD')
    
            overall_meanzbins.append(meanzinterp[0]) # no longer in 5 bins
            # we then multiply the sensitivity from Jeff's code by the std of the i-band mag
            #to get the std of the z 
            
            stdz = [float(dzdminterp)*float(bd[list(bd.keys())[0]].summary_values['Rms']) ] 
            # using chain rule and squaring deriv
 
            # Hard coding this to check interpolation
            #stdz = [float(0.03) ]
            
            overall_stdzbins.append(stdz[0])
            clbias, meanz_use = compute_Clbias(meanzinterp,stdz)
            overall_clbias.append(clbias[0])
            meanz_usecl.append(meanz_use[0])
            
            
    df = pd.DataFrame(list(zip(overall_names, overall_years, overall_meds, overall_means, overall_std, overall_iqr, overall_meanzbins,overall_stdzbins, overall_clbias, meanz_usecl)), 
                  columns=['Strategy', 'Year', 'Median i-band depth', 'Mean i-band depth', 'Std i-band depth', 'IQR i-band depth', 'Mean z', 'Std z','Clbias','Used meanz'])
    return df

def compute_dzfromdm(zbins, imag, year, dzname):


    if dzname='SRD':
        #dzname = Jeff's implementation
            
            print('results_%i.feather'%year)
            zgrid = pd.read_feather('results_%i.feather'%year)  
            dz = zgrid['meanz']-zgrid['true_meanz']
            dz = dz.values
            meanz = zgrid['meanz'].values # no longer in 5 bins
            m5s = zgrid['m5'].values
            # calculate the sample i magnitude limit used to generate this file
            ilim_year = 25.3+1.25*np.log10(year/10) 
            # calculate the derivative of mean z with respect to m5    
            dzdm = np.gradient(meanz)/np.gradient(m5s)
            # assuming here that imag is the actual depth of the imaging
            #       (so it should be compared to what is tabulated in m5s)
            dzdminterp = np.interp(imag, m5s, dzdm)
            dzdminterp = 0.0026 # taken from Jeff's average 
            meanzinterp = [np.interp(imag, m5s, meanz)]
    return dzdminterp, meanzinterp

# Define combined plotting routine - base it on Renee's
def combined_metric_plots(use_run_name_vec, use_opsim_fname_vec, 
                          use_metric=maf.ExgalM5(), year=10, use_color_min=None, use_color_max=None,nside=64):
    # use_run_name_vec says which OpSim DBs we want to use - will also be used for labels
    # use_opsim_fname_vec says where they live, e.g. one item might be `/global/cfs/cdirs/lsst/groups/CO/rubin_sim/sim_baseline/baseline_v2.1_10yrs.db`
    if use_color_min is not None and use_color_max is not None:
        plot_dict={"color_min": use_color_min, "color_max": use_color_max, "x_min": use_color_min, "x_max": use_color_max}
    else:
        plot_dict=None
    days = year*365.3
    constraint_str='filter="i" and note not like "DD%" and night <= XX and note not like "twilight_neo" '
    constraint_str = constraint_str.replace('XX','%d'%days)
    print(constraint_str)
    
    bg_list = []
    bd_list = []
    overall_plot_dict = {}
    color_list = ["k", "r", "b", "c", "g", "o", "m", "y"] # hopefully long enough to handle everything
    for i in range(len(use_run_name_vec)):
        use_run_name = use_run_name_vec[i]
        use_opsim_fname = use_opsim_fname_vec[i]
        print(use_run_name, use_opsim_fname)
        depth_map_bundle = maf.MetricBundle(
            metric=use_metric,
            slicer=maf.HealpixSubsetSlicer(nside=nside, use_cache=False, hpid=np.where(map_labels == "lowdust")[0]),
            constraint=constraint_str,
            run_name=use_run_name,
            summary_metrics=[maf.MedianMetric(), maf.MeanMetric(), maf.RmsMetric()],
            plot_dict=plot_dict
        )
        #print('Bundle diagnostics',depth_map_bundle.run_name, depth_map_bundle.metric.name, 
        #      depth_map_bundle.info_label, depth_map_bundle.slicer.slicer_name, depth_map_bundle.file_root)
        
        bd = maf.metricBundles.make_bundles_dict_from_list([depth_map_bundle])
        print(bd[list(bd.keys())[0]].summary_values)
        #bgroup = maf.MetricBundleGroup(
        #    bd, use_run_name, out_dir="./"
        #)
        #bgroup.run_all()
        
        #bg_list.append(bgroup)
        #bd_list.append(bd)
        #overall_plot_dict[use_run_name] = color_list[i]
    
    #ph = maf.PlotHandler()
    #ph.set_metric_bundles(bg_list)
    #ph.plot(plot_dicts=overall_plot_dict)
    #for i in range(len(bd_list)): print(use_opsim_fname_vec[i], year, bd_list[i][list(bd_list[i].keys())[0]].summary_values)
    #for tmpbd in bd_list: print(tmpbd[list(tmpbd.keys())[0]].summary_values)
 
def use_zbins(meanz_vals, figure_9_mean_z=np.array([0.2, 0.4, 0.7, 1.0]),  figure_9_width=0.2):
    max_z_use = np.max(figure_9_mean_z)+2*figure_9_width
    use_bins = meanz_vals < max_z_use
    #print(use_bins, 'use_bins')
    return use_bins

def compute_Clbias(meanz_vals,scatter_mean_z_values,figure_9_mean_z=np.array([0.2, 0.4, 0.7, 1.0]), figure_9_Clbias =np.array([1e-3, 2e-3, 5e-3, 1.1e-2]),figure_9_width=0.2,figure_9_mean_z_scatter = 0.02):
    import numpy as np
    mzvals= np.array([float(mz) for mz in meanz_vals])
    sctz = np.array([float(sz)for sz in scatter_mean_z_values])
    
    fit_res = np.polyfit(figure_9_mean_z, figure_9_Clbias, 2)
    poly_fit = np.poly1d(fit_res)
    use_bins = use_zbins(meanz_vals,figure_9_mean_z, figure_9_width)

    mean_z_values_use = mzvals[use_bins]
    sctz_use = sctz[use_bins]

    Clbias = poly_fit(mean_z_values_use)
    rescale_fac =  sctz_use / figure_9_mean_z_scatter
    Clbias *= rescale_fac
    fit_res_bias = np.polyfit(mean_z_values_use, Clbias, 1)
    poly_fit_bias = np.poly1d(fit_res_bias)

    return poly_fit_bias(mean_z_values_use), mean_z_values_use


# NOT USED BELOW, READ DIRECTLY FROM JEFF'S FILES
# def n_of_i_func(imin=17,imax=28,ni=101, zmin=0,zmax=4,nz=401,n_mc=300000):

#     # set up to make interpolation tables
#     # can change the i magnitude range of objects or the z
#     # range to consider in this cell
#     deltai = (imax-imin)/(ni-1)
#     deltaz = (zmax-zmin)/(nz-1)
#     ival=np.linspace(imin,imax,ni)
#     zvals=np.linspace(zmin,zmax,nz)

#     # generate a normalized cumulative distribution 
#     # corresponding to n(<i) from the SRD
#     # for the Monte Carlo we want the values to run from 0 to 1

#     n_of_i = 42.9*0.88*10**(0.359*(ival-25))
#     n_of_i = n_of_i / n_of_i[ni-1]
 
#     return ival, n_of_i,zvals


# def generate_zdistribution(n_mc=3000000,imin=17,imax=28,ni=101, zmin=0,zmax=4,nz=401,filename='zdist.pkl'):
#     ''' Code from Jeff Newman to generate the distributions of objects with given limiting magnitude. 
#     Jeff notes that you need > 1M MC iterations for good errorbars - 
#     but we are reducing the default value here for speed.'''
#     import pandas as pd
    
#     ival,n_of_i,zvals = n_of_i_func(imin,imax,ni,zmin,zmax,nz,n_mc)
#     # generate random true i magnitudes up to imax following the 
#     # distribution calculated in the above cell
#     true_i = np.interp(np.random.random_sample(n_mc),n_of_i,ival)

# # array to contain the true redshift for each MC object
#     true_z = np.zeros_like(true_i)

#     #based on the drawn i magnitude for each object,
#     #   draw a random z from p(z,i)
#     # i was hoping to avoid calculating all this for every
#     # object but failed to get interpolation to work.  This 
#     # does ok but is a bit slower than i'd like

#     for idx,imag in enumerate(true_i):
#         z0 = 0.246 + 0.025*(imag-24.1)
#         dndi = (0.868296*10**(0.359*(imag-25))
#          *np.exp(-(zvals/z0)**0.92)*zvals**3)/( (z0**2)*(zvals/z0)**0.08) 
#         dndi += 31.2069*10**(0.359*(imag-25))*np.exp(-(zvals/z0)**0.92)*zvals**2
#         dndi[np.isfinite(dndi) == 0] = 0.
    
#         cdndi = np.cumsum(dndi)
#         cdndi = cdndi / cdndi[nz-1]
#         true_z[idx] = np.interp(np.random.random_sample(1),
#           cdndi,zvals)
        
#     d= {'i':true_i,'zarr':true_z}

#     catalog = pd.DataFrame(data=d)
#     catalog.to_pickle(filename)
#     return catalog
    

# def compute_deltaz(generate_zdist=False,zdistfile='zdist.pkl',n_mc=30000,imin=17,imax=28,ni=101,zmin=0,zmax=4,nz=401,imag=25.3,catalog_mc=100000,flux_var=0.01,m5=26):
#     from photerr import LsstErrorModel

#     if generate_zdist:
#         catalog = generate_zdistribution(n_mc,imin,imax,ni,zmin,zmax,nz,filename=zdistfile)
#     else:
#        catalog = pd.read_pickle(zdistfile)

#     catalog_rep = catalog.copy()
#     errModel = LsstErrorModel(nYrObs=1,nVisYr={'i':1},m5={'i':float(m5)})

#     tmpcatalog = errModel(catalog_rep, random_state=np.random.randint(1,catalog_mc))
#     fluxes = 10**(-0.4*(tmpcatalog['i'] - 27))
#     fluxerrs = tmpcatalog['i_err']*np.log(10)/2.5*fluxes

#     noisy_i_flux = fluxes + fluxerrs*np.random.normal(size=tmpcatalog.count()[0])
#     # don't let fluxes go negative: limit corresponds to magnitude = 32
#     noisy_i_flux = np.maximum(noisy_i_flux, flux_var)
#     noisy_i = 27 - 2.5*np.log10(noisy_i_flux)

#     meanz_imag = np.mean(tmpcatalog[noisy_i < imag]['zarr'])
#     number_imag = np.sum(noisy_i < imag)
#     meanz_imag_error = np.std(tmpcatalog[noisy_i < imag]['zarr'])/np.sqrt(number_imag)
#     number_imag_error = np.sqrt(number_imag)

#     return meanz_imag, number_imag, meanz_imag_error, number_imag_error


# def grid_deltaz(num_m5s=26,m5min=28.25,m5max=25.75,imag=25.3,catalog_mc=100000,flux_var=0.01, 
#                 generate_zdist=False,zdistfile='results_bright.pkl',
#                 n_mc=3000000,imin=17,imax=28,ni=101,zmin=0,zmax=4,nz=401):

#     if generate_zdist:
#         catalog = generate_zdistribution(n_mc,imin,imax,ni,zmin,zmax,nz,filename=zdistfile)
#     else:
#        catalog = pd.read_pickle(zdistfile)


   
#     m5s = np.linspace(m5max,m5min,num_m5s)
#     meanz_imag = np.zeros_like(m5s)
#     meanz_imag_error = np.zeros_like(m5s)
    

#     number_imag = np.zeros_like(m5s)
#     number_imag_error = np.zeros_like(m5s)

#     true_z = catalog['zarr']
#     true_i = catalog['i']

#     for idx,m5 in enumerate(m5s):
#         meanz_imag[idx], number_imag[idx], meanz_imag_error[idx], number_imag_error[idx] = compute_deltaz(generate_zdist=False,zdistfile='zdist.pkl',
#                                                                                                           n_mc=n_mc,imin=imin,imax=imax,ni=ni,zmin=zmin,
#                                                                                                           zmax=zmax,nz=nz,imag=imag,catalog_mc=catalog_mc,
#                                                                                                           flux_var=flux_var,m5=m5)
        
#     meanz_imag_true = np.mean(catalog[true_i < imag]['zarr'])
#     number_imag_true = np.sum(true_i < imag)

#     d= {'true_meanz':meanz_imag_true,'true_n':number_imag_true,'m5s':m5s,'meanz':meanz_imag,
#      'number': number_imag, 'meanz_err':meanz_imag_error,'number_err':number_imag_error }
#     outputs = pd.DataFrame(data=d)

#     outputs.to_pickle('deltazgrid_outputs_dataframe.pkl')
#     return outputs
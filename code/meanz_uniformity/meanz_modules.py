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

def get_3x2pt_metric(use_run_name, use_opsim_fname, year=10, nside=64, return_map=False, my_hpid=None):
    # use_run_name says which OpSim DB we want to use, e.g. `baseline_v2.1_10yrs` - will also be used for labels
    # use_opsim_fname says where it lives, e.g. `/global/cfs/cdirs/lsst/groups/CO/rubin_sim/sim_baseline/baseline_v2.1_10yrs.db`
    surveyAreas = SkyAreaGenerator(nside=nside)
    map_footprints, map_labels = surveyAreas.return_maps()
    days = year*365.25
    use_filter = "i"
    # Here the constraint on use of i-band data, exclusion of DDFs, time limitations, and avoiding twilight exposures 
    constraint_str='filter="YY" and note not like "DD%" and night <= XX and note not like "twilight_near_sun" '
    constraint_str = constraint_str.replace('XX','%d'%days)
    constraint_str = constraint_str.replace('YY','%s'%use_filter)
    
    ThreebyTwoSummary = maf.StaticProbesFoMEmulatorMetric(nside=nside, metric_name="3x2ptFoM")

    # Decide what summary statistics we want
    my_summary_stats = [maf.MedianMetric(), maf.MeanMetric(), maf.RmsMetric(), maf.PercentileMetric(percentile=25), maf.PercentileMetric(percentile=75), ThreebyTwoSummary]
    
    # First, define a MetricBundle object.
    if my_hpid is None: 
        use_slicer = maf.HealpixSubsetSlicer(nside=nside, use_cache=False, hpid=np.where(map_labels == "lowdust")[0])
    else:
        use_slicer = maf.HealpixSubsetSlicer(nside=nside, use_cache=False,
                                             hpid=np.intersect1d(np.where(map_labels == "lowdust")[0], my_hpid))
    depth_map_bundle = maf.MetricBundle(
        metric=maf.ExgalM5(),
        # Exclude the galactic plane
        slicer=use_slicer,
        constraint=constraint_str,
        run_name=use_run_name,
        summary_metrics=my_summary_stats
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

def my_total_power_metric(map, ell_max=30, return_functions=False):
    use_map = map.copy()
    mean_val = np.mean(map[map>0])
    use_map[map>0] -= mean_val
        
    cl = hp.anafast(use_map)
    ell = np.arange(np.size(cl))
    total_power = np.sum((2*ell[ell<ell_max]+1)*cl[ell<ell_max])
    if return_functions:
        return total_power, ell, cl
    else:
        return total_power

# A utility to plot summary stats for strategies as a function of year, given a dataframe from the above routines.
def plot_metric_by_year(df, stat_name,years=None,filter='i', y_axis_label=None,ylog=False, compare_to=None):

    if years!=None:
        year_vals = years
    else:
        year_vals = np.array(list(set(df['Year'])))

    cols = ['lightcoral', 'mediumpurple','deepskyblue','teal','forestgreen','r','g','k','c','m']

    # print(year_vals, 'years')
    # filter_use=filter
    # print(filter_use, 'filter')
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
    for sc,s in enumerate(strategies):
        metricvals=[]
        comparevals=[]
        #yvals = np.array([float(val) for val in df[stat_name][df['Strategy']==s]])
        print(s)
        for i,year in enumerate(year_vals):
            if filter=='combined':
                tmpval=df[stat_name][df['Strategy']==s][df['Year']==year].values[0]
                metricvals.append(tmpval)
                if compare_to!=None:
                    tmpval = df[stat_name][df['Strategy']==compare_to][df['Year']==year].values[0]
                    comparevals.append(tmpval)
            else:
                tmpval=df[stat_name][df['Strategy']==s][df['Year']==year].values[0][filter]
                metricvals.append(tmpval)
                if compare_to!=None:
                    tmpval = df[stat_name][df['Strategy']==compare_to][df['Year']==year].values[0][filter]
                    comparevals.append(tmpval)
            

        
        metricvals=np.array(metricvals)
        if compare_to!=None:
            ax.plot(year_vals+1+offsets[offset_index], (metricvals-comparevals)/comparevals, color=cols[sc],label=s) 
            plt.ylabel(f'Fractional change in {y_axis_label}')
        else:   
            ax.plot(year_vals+1+offsets[offset_index], metricvals, color=cols[sc],label=s)
            plt.ylabel(y_axis_label)
        offset_index += 1
    plt.xlabel('Year')
    
    if ylog: plt.yscale('log')
    plt.legend()
    plt.show()
    
# A utility to plot summary stats for strategies as a function of year, given a dataframe from the above routines.
def plot_meanz_metrics_by_year(df, years, filter='i',num_bins=5,y_axis_label=None):
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
        if (sy==0) and (filter=='combined'):
            axs[0][sy].set_title(f'Combined offsets: Year {year+1}', fontsize=15)
        elif (sy==0):
            axs[0][sy].set_title(f'Offsets for filter {filter}: Year {year+1}', fontsize=15)
        else:
            axs[0][sy].set_title(f'Year {year+1}', fontsize=15)
        #axs[1][sy].set_xlabel('Tomographic Bin', fontsize=15)
        axs[2][sy].set_xlabel('Tomographic Bin', fontsize=15)
        bins=np.arange(num_bins)
        for scount,s in enumerate(strategies):

            if filter=='combined':
                meanz=df['Average mean z'][df['Strategy']==s][df['Year']==year].values[0]
                stdz=df['Combined std z'][df['Strategy']==s][df['Year']==year].values[0]
            else: 
                meanz=df['Mean z'][df['Strategy']==s][df['Year']==year].values[0][filter]
                stdz=df['Std z'][df['Strategy']==s][df['Year']==year].values[0][filter]
                
            if (sy==len(year_vals)-1):
                axs[0][sy].plot(bins+offset*scount,meanz,label=s,color=cols[scount])
                axs[1][sy].plot(bins+offset*scount,stdz, label=s,color=cols[scount])
            else:
                axs[0][sy].plot(bins+offset*scount,meanz, color=cols[scount])
                axs[1][sy].plot(bins+offset*scount,stdz, color=cols[scount])
        
            if filter=='combined':
                nbins_use = np.shape(df['Combined Mean z use'][df['Strategy']==s][df['Year']==year].values[0])[0]
                clbias = df['Combined Clbias'][df['Strategy']==s][df['Year']==year].values[0]
            else:
                nbins_use = np.shape(df['Used meanz'][df['Strategy']==s][df['Year']==year].values[0][filter])[0]
                clbias = df['Clbias'][df['Strategy']==s][df['Year']==year].values[0][filter]

            bins_use=np.arange(nbins_use)
            if (sy==len(year_vals)-1):
                axs[2][sy].plot(bins_use+offset*scount,clbias, label=s,color=cols[scount])
            else:
                axs[2][sy].plot(bins_use+offset*scount,clbias, color=cols[scount])

    #     
        #axs[sy].set_yscale('log')


    axs[0][0].set_ylabel('Mean z', fontsize=15)
    axs[1][0].set_ylabel('Std z', fontsize=15)
    axs[2][0].set_ylabel('Cl bias', fontsize=15)
    axs[1][sy].legend(loc='upper left',fontsize=10)
    #axs[0].show()
# A utility to plot summary stats for strategies as a function of year, given a dataframe from the above routines.


def get_year_by_year_metrics(year_list, name_list, sim_list):

    metricList=[]
    filter_list=["u","g","r","i"]#,"z","y"]
    zbins=5
    filts = dict(zip(filter_list, [None]*len(filter_list)))
    filtsz = dict(zip(filter_list, [[None]*zbins]*len(filter_list)))
    keyList = ['Strategy','Year','Mean depth','Median depth','Std depth','IQR depth','Mean z','Std z','Clbias','Used meanz']
    metricDict = {} #dict(zip(keyList, [None]*len(keyList)))
    for s in ['Mean depth','Median depth','Std depth','IQR depth']:
        metricDict[s] = dict(zip(filter_list, [None]*len(filter_list)))
    for s in ['Mean z','Std z','Clbias','Used meanz']:
        metricDict[s]=dict(zip(filter_list, [[None]*zbins]*len(filter_list)))


    for year in year_list:
        
        for i,sim in enumerate(sim_list):
            tmpmetricDict = {} #dict(zip(keyList, [None]*len(keyList)))

            

            for s in ['Mean depth','Median depth','Std depth','IQR depth']:
                tmpmetricDict[s] = dict(zip(filter_list, [None]*len(filter_list)))

            for s in ['Mean z','Std z','Clbias','Used meanz']:
                tmpmetricDict[s]=dict(zip(filter_list, [[None]*zbins]*len(filter_list)))

            totdz=np.zeros(zbins)
            avmeanz=np.zeros(zbins)
            combinedclbias=np.zeros(zbins)
            for filter_ind,use_filter in enumerate(filter_list):

                tmpmetricDict['Strategy']=name_list[i] # strategy name
                tmpmetricDict['Year']=year
                bgroup, bd = metric_plots(name_list[i], sim_list[i], year=year, use_filter=use_filter)
                tmpmetricDict['Mean depth'][use_filter]=bd[list(bd.keys())[0]].summary_values['Mean'] 
                tmpmetricDict['Median depth'][use_filter]=bd[list(bd.keys())[0]].summary_values['Median']
                tmpmetricDict['Std depth'][use_filter]=bd[list(bd.keys())[0]].summary_values['Rms'] # rms of the mags
                tmpmetricDict['IQR depth'][use_filter]= bd[list(bd.keys())[0]].summary_values['75th%ile']-bd[list(bd.keys())[0]].summary_values['25th%ile'] 
                dzdminterp, meanzinterp=compute_dzfromdm(zbins, filter_ind,year, 'JQ')
                tmpmetricDict['Mean z'][use_filter] = meanzinterp 
                # we then multiply the sensitivity from Jeff's code by the std of the mag
                # to get the std of the z 
                stdz = [float(np.abs(dz))*float(bd[list(bd.keys())[0]].summary_values['Rms']) for dz in dzdminterp]
                tmpmetricDict['Std z'][use_filter] = stdz 
                clbias, meanz_use = compute_Clbias(meanzinterp,stdz)
                tmpmetricDict['Clbias'][use_filter]=clbias
                tmpmetricDict['Used meanz'][use_filter]=meanz_use

                totdz+=stdz
                avmeanz+=meanzinterp
           
            #print(avmeanz)
            tmpmetricDict['Combined std z'] = totdz
            tmpmetricDict['Average mean z']= avmeanz/len(filter_list)
            tmpclbias,tmpmeanz_use=compute_Clbias(avmeanz/len(filter_list),totdz)
            tmpmetricDict['Combined Clbias']=tmpclbias
            tmpmetricDict['Combined Mean z use']=tmpmeanz_use
            y1rat,y10rat = compare_Clbias_DESCreg(tmpclbias)

            tmpmetricDict['Y1 ratio']=y1rat
            tmpmetricDict['Y10 ratio']=y10rat

            metricList.append(tmpmetricDict)

    df = pd.DataFrame(metricList)    
    return df

def compute_dzfromdm(zbins, band_ind, year, dzname):

    if dzname=='JQ':
        #dzname = Jeff's implementation
        deriv = pd.read_pickle('uniformity_pkl/meanzderiv.pkl')
        zvals = pd.read_pickle('uniformity_pkl/meanzsy%i.pkl'%(year+1))
        meanz = zvals[0:zbins,band_ind,5] # need to think through this, 
        #currently taking an arbitrary delta index from Jeff's code
        dz = deriv[year,band_ind,0:zbins]
        dzdminterp = np.abs(dz)
        meanzinterp = meanz
        #print(meanzinterp,  'meanz, checking implementation')
        #print(dzdminterp,'dzdm, checking implementation')
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
    # print(constraint_str)
    # print('----------------')
    
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
    
    # print(mzvals, 'mzvals in computclbias')
    # print(sctz, 'sctz in computclbias')
    fit_res = np.polyfit(figure_9_mean_z, figure_9_Clbias, 2)
    poly_fit = np.poly1d(fit_res)
    use_bins = use_zbins(meanz_vals,figure_9_mean_z, figure_9_width)
    # print('use bins ', use_bins)
    mean_z_values_use = mzvals[use_bins]
    sctz_use = sctz[use_bins]

    Clbias = poly_fit(mean_z_values_use)
    rescale_fac =  sctz_use / figure_9_mean_z_scatter
    Clbias *= rescale_fac
    fit_res_bias = np.polyfit(mean_z_values_use, Clbias, 1)
    poly_fit_bias = np.poly1d(fit_res_bias)

    return poly_fit_bias(mean_z_values_use), mean_z_values_use


def compare_Clbias_DESCreg(clbias):

    y10_req = 0.003
    y1_goal = 0.013

    clbiastot = np.max(clbias)
    y10ratio = clbiastot/y10_req
    y1ratio = clbiastot/y1_goal

    return(y1ratio,y10ratio)

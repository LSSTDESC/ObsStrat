import numpy as np
import pandas as pd

ylabel = '(metric - baseline)/baseline'
# baseline = 'baseline_v1.4_10yrs.db'
baseline = 'baseline_v1.5_10yrs'
# baseline = 'baseline_nexp1_v1.6_10yrs.db'

metrics_to_drop = ['kn_GW170817_error', 
                    'kn_Pop_kNe_error',
                    'sl_td_all_bands_distance_precision',
                    'sl_std_num_lensed_SNe_Ia_in_specific_galaxy_clusters ',
                    'pz_y10_zlow_bias',
                    'pz_y10_zhigh_bias']

metrics_to_drop += [
    'static_median Y3 $i$-band coadded depth in effective survey area',
    'static_median Y6 $i$-band coadded depth in effective survey area',
    'static_stddev of Y3 $i$-band coadded depth in effective survey area',
    'static_stddev of Y6 $i$-band coadded depth in effective survey area',
    'static_Y3 effective survey area',
    'static_Y6 effective survey area'
] # Removed at the request of Humna and Eric

def read_metrics(flname, baseline='baseline_v1.5_10yrs', dropna=True,
                 add_samefilt=False, drop_baseline=False):
    metrics = pd.read_csv(flname, index_col=0)
    metric_display_df = pd.read_csv('metric_display_names.csv', index_col=0)
    

    if add_samefilt:
        run_to_add = 'baseline_samefilt_v1.5_10yrs.db'
        metrics_outstanding = pd.read_csv('metrics_outstanding.csv', index_col=0)
        if run_to_add not in metrics.index:
            this_df = metrics_outstanding.loc[[run_to_add]]
            metrics = pd.concat((metrics, this_df))

    

    runs_to_drop = ['testrolling', 'stuck_rolling']
    runs_to_drop += ['alt_base_v1.6_10yrs', 'rolling_fpo_v1.6_10yrs'] # v1.6

    if 'Problematic' in metrics.columns:
        # metrics.loc['u60_v1.5_10yrs', 'Problematic'] = pd.NA
        obsolete_runs = metrics.index[~metrics.Problematic.isna().values]
        runs_to_drop += list(obsolete_runs)
        metrics = metrics.drop('Problematic', axis=1)
    
    if 'group' in metrics.columns:
        metrics = metrics.drop('group', axis=1)

    convert_to_flux_cols = ['static_median Y10 $i$-band coadded depth in effective survey area',
                            'static_stddev of Y10 $i$-band coadded depth in effective survey area',
                            'static_median Y1 $i$-band coadded depth in effective survey area',
                            'static_stddev of Y1 $i$-band coadded depth in effective survey area']


    one_minus_cols = ['pz_y10_zlow_fout', 'pz_y10_zhigh_fout']
                            
    convert_to_precision_cols = ['static_median Y10 $i$-band coadded depth in effective survey area',
                                'static_stddev of Y10 $i$-band coadded depth in effective survey area',
                                'static_median Y1 $i$-band coadded depth in effective survey area',
                                'static_stddev of Y1 $i$-band coadded depth in effective survey area',
                                'pz_y10_zlow_stdd',
                                'pz_y10_zhigh_stdd',
                                'pz_y10_zhigh_bias', 
                                'pz_y10_zlow_bias'
                                ]

    # sqrt_cols = ['sn_normal_nsn_tot']
    sqrt_cols = []

    # new_display_names = {
    # 'static_median Y10 $i$-band coadded depth in effective survey area': 'Median depth',
    # 'static_stddev of Y10 $i$-band coadded depth in effective survey area': 'Inverse variance of depth',
    # 'static_Y10 effective survey area': 'Effective survey area',
    # 'lss_diagnostic FOM Y1': 'LSS Systematics FoM Y1',
    # 'lss_diagnostic FOM Y10': 'LSS Systematics FoM Y10',
    # 'lss_Y1 Ngal (0.66<z<1.0)': 'Ngal Y1',
    # 'lss_Y10 Ngal (0.66<z<1.0)': 'Ngal Y10',
    # 'sl_num_lensed_SNe_Ia_with_good_time_delays_lensed_by_galaxies': 'No. SNIa lensed by galaxies',
    # 'sl_num_lensed_quasars_with_good_time_delays': 'No. lensed quasars',
    # 'sl_num_lensed_SNe_Ia_in_specific_galaxy_clusters': 'No. SNIa lensed by clusters',
    # 'wl_systematics': 'WL systematics proxy',
    # 'sn_faint_nsn_tot': 'No. faint SNIa', 
    # 'sn_faint_zmax': 'z_max of faint SNIa', 
    # 'sn_normal_nsn_tot': 'No. normal SNIa',
    # 'sn_faint_snr_rate': 'Faint SNIa SNR rate', 
    # 'sn_faint_zlim_r': 'Faint SNIa z_lim', 
    # 'sn_normal_zmax': 'z_max of normal SNIa',
    # 'pz_y10_zlow_fout': 'PZ inlier fraction low z',
    # 'pz_y10_zlow_stdd': 'PZ inverse variance low z',
    # 'pz_y10_zlow_bias':	'PZ inverse bias low z',
    # 'pz_y10_zhigh_fout': 'PZ inlier fraction high z',	
    # 'pz_y10_zhigh_stdd': 'PZ inverse variance high z',	
    # 'pz_y10_zhigh_bias': 'PZ inverse bias high z'
    # }

    # words_to_replace = {
    #                     'outlier': 'inlier',
    #                     'stddev': 'precision',
    #                     'bias': 'inverse bias$^2$'
    # }
    words_to_replace = {}

    drop_intervening_years = True

    
    metric_display_df = metric_display_df.drop(metrics_to_drop, axis=1)
    metrics_display_name = metric_display_df.loc['DISPLAY_NAME']

    # for k in new_display_names.keys():
    #     metrics_display_name.loc[k] = new_display_names[k]
    for metric in metrics_display_name.index:
        s = metrics_display_name.loc[metric]
        for k in list(words_to_replace.keys()):
            if k in s:
                new_s = s.replace(k, words_to_replace[k])
                metrics_display_name.loc[metric] = new_s


    # metrics = metrics.drop(['DISPLAY_NAME'])
    metrics = metrics.drop(metrics_to_drop, axis=1)
    if 'group' in metrics.columns:
        metrics.drop('group', axis=1)

    for c in metrics.columns:
        metrics[c] = metrics[c].astype('float')
        
        if drop_intervening_years:
            if 'static' in c and 'Y10' not in c and 'Y1' not in c:
                metrics = metrics.drop(c, axis=1)
        
        # if c=='wl_fom':
        #     metrics = metrics.rename(columns={'wl_fom':'detf_fom'})
        #     metrics_display_name = metrics_display_name.rename({'wl_fom':'detf_fom'})

    for i in metrics.index:
        for r in runs_to_drop:
            if r in i:
                if i in metrics.index: # may already be dropped
                    metrics = metrics.drop(i, axis=0)
    if dropna:
        metrics = metrics.dropna(axis=1)

    # Do this before all the conversions?
    baseline_values = metrics.loc[baseline].to_dict()

    for c in convert_to_flux_cols:
        if c in metrics.columns:
            metrics[c] = 10**(-metrics[c]/2.5)

    for c in one_minus_cols:
        if c in metrics.columns:
            metrics[c] = 1 - 2*metrics[c]

    for c in sqrt_cols:
        if c in metrics.columns:
            metrics[c] = np.sqrt(metrics[c])

    
    if 'inter_night_gap_all_filters' in metrics.columns:
        metrics['median_num_obs_per_night'] = 1/metrics['inter_night_gap_all_filters']
        metrics.drop('inter_night_gap_all_filters', axis=1)

    for c in convert_to_precision_cols:
        if c in metrics.columns:
            metrics[c] = 1/metrics[c]**2

    metrics_rel = metrics.copy()

    # baseline_values = metrics.loc[baseline].to_dict()

    for c in metrics_rel.columns:
        metrics_rel[c] = (metrics_rel[c] - metrics_rel.loc[baseline,c])/metrics_rel.loc[baseline,c]
        # metrics_rel[c] = (metrics_rel[c])/metrics_rel.loc[baseline,c]
    
    if drop_baseline:
        metrics_rel = metrics_rel.drop(baseline)
        
    return metrics_rel, metrics_display_name, baseline_values

def get_metrics_by_probe(metrics_rel):
    probes_metrics = {}
    prefixes = ['static', 'lss', 'wl', 'sl', 'sn', 'kn', 'pz']

    for pre in prefixes:
        cols = [c for c in metrics_rel.columns if pre+'_' in c]
        #cols = ['group'] + cols
        probes_metrics[pre] = metrics_rel[cols]
    
    return probes_metrics
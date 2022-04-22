import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import george
from scipy.optimize import minimize
import seaborn as sns
import read_metrics

which_probes = 'all'

if which_probes == 'all':
    probes_translations = {'fom':'wl_fom',
                        'lss':'lss_diagnostic FOM Y10',
                        'wl':'wl_systematics',
                        'sn':'sn_normal_nsn_tot',
                        'snz': 'sn_normal_zmax',
                        'sl':'sl_num_lensed_SNe_Ia_with_good_time_delays_lensed_by_galaxies',
                        'kn':'kn_Pop_kNe_counts',
                        'pz': 'pz_y10_zlow_stdd'
                        }

    colors = {'fom':'C0',
            'lss':'C1',
            'wl':'C2',
            'sn':'C3',
            'sl':'C4',
            'snz':'C5',
            'kn':'C6',
            'pz':'C9'
        }

    markers = {'fom':'o',
            'lss':'^',
            'wl':'s',
            'sn':'*',
            'snz':'P',
            'sl':'D',
            'kn':'X',
            'pz':'+'
        }

elif which_probes == 'pz':

    probes_translations = {'pz_y10_zlow_fout': 'pz_y10_zlow_fout',
                            'pz_y10_zlow_stdd': 'pz_y10_zlow_stdd',
                            'pz_y10_zlow_bias': 'pz_y10_zlow_bias',
                            'pz_y10_zhigh_fout': 'pz_y10_zhigh_fout',
                            'pz_y10_zhigh_stdd': 'pz_y10_zhigh_stdd',
                            'pz_y10_zhigh_bias': 'pz_y10_zhigh_bias'
                            }

    colors = {'pz_y10_zlow_fout':'C0',
                'pz_y10_zlow_stdd':'C1',
                'pz_y10_zlow_bias':'C2',
                'pz_y10_zhigh_fout':'C3',
                'pz_y10_zhigh_stdd':'C4',
                'pz_y10_zhigh_bias':'C6'
            }

    markers = {'pz_y10_zlow_fout':'o',
                'pz_y10_zlow_stdd':'^',
                'pz_y10_zlow_bias':'s',
                'pz_y10_zhigh_fout':'*',
                'pz_y10_zhigh_stdd':'D',
                'pz_y10_zhigh_bias':'X'
            }

t0 = {'lw':2, 'alpha':1}
t1 = {'lw':1, 'alpha':0.8}

def make_bar_graph(this_df, metrics_display_name, orientation='vertical', 
                    ylabel='', fname='', 
                    descriptions=[], plot_data=[],
                    baseline_values={},
                    legend_location='best', plot_type='presentation',
                    figsize=(10,15), expand_ylim_by=[1.0, 1.0],
                    expand_xlim_by=[1.0,1.0],
                    include_baseline_values=False,
                    include_colour_bands=False,
                    include_descriptions=False,
                    use_short_names=True,
                    shaded_pairs=False):
    plt.close('all')
    if 'group' in this_df.columns or 'Problematic' in this_df.columns:
        this_df = this_df[this_df.columns[1:]]


    #cols = plt.rcParams['axes.prop_cycle'].by_key()['color']
    cols = [(1,1,1), (0.6,0.6,0.6)] # This only allows two colours now

    # if orientation == 'vertical':
    #     figsize = (10,15)
    # else: figsize = (10, 8)

    fig = plt.figure(figsize=figsize)
    wd = 0.8/(len(this_df.columns))
    x = np.arange(len(this_df.index))
    x1extra = np.arange(len(this_df.index)+1)
    

    maxi = this_df.max().max()
    mini = this_df.min().min()

    lim = max(abs(mini), maxi)
    lims = [-lim*expand_xlim_by[0], lim*expand_xlim_by[1]]
    index_lims = [-0.2, len(this_df.index)+0.5]

    if orientation == 'vertical':
        
        for i in np.arange(len(this_df.index)+1):
            plt.plot(lims, [i-0.1, i-0.1], '--k', alpha=0.3)
    else:
        plt.plot([0, len(this_df.index)+1], [0,0], 'k', alpha=0.3)

    if len(plot_data) != 0:
        for i in range(len(x)):
            this_sim = this_df.index[i]
            if this_sim in plot_data.index:
                color = cols[plot_data.loc[this_sim, 'color']]
                if include_colour_bands:
                    plt.axhspan(x1extra[i]-0.1,x1extra[i+1]-0.1, 
                                    facecolor=color, alpha=0.1)
                
                annotation = plot_data.loc[this_sim, 'annotation']
                annotation = annotation.replace('\\n', '\n')
                xoffset = plot_data.loc[this_sim, 'xoffset']
                yoffset = plot_data.loc[this_sim, 'yoffset']
                plt.text(xoffset, i + yoffset, annotation)
            
    bar_colours = ['C%d' %i for i in range(10)]
    k=0
    for c in this_df.columns:
        lbl = metrics_display_name[c]
        if plot_type == 'paper':
            lbl = lbl.replace('<', '$\\textgreater$')
            lbl = lbl.replace('>', '$\\textless$')
            lbl = lbl.replace('_{', '$_{')
            lbl = lbl.replace('}', '}$')
        if len(list(baseline_values.keys())) != 0 and include_baseline_values:
            lbl += ' (%.3g)' %baseline_values[c]

        if shaded_pairs:
            bar_colour = bar_colours[k//2]
            if k%2 != 0:
                hatch = '////'
            else:
                hatch = None
        else:
            bar_colour = bar_colours[k]
            hatch = None
        edgecolor = [1]*3
        if orientation == 'vertical':
            plt.barh(x+k*wd, this_df[c], height=wd, align='edge', label=lbl,
                     facecolor=bar_colour, hatch=hatch, alpha=0.99,
                     edgecolor=edgecolor, lw=0)
        else:
            plt.bar(x+k*wd, this_df[c], width=wd, align='edge', label=lbl,
                     facecolor=bar_colour, hatch=hatch, alpha=0.99,
                     edgecolor='k')

        k+=1
    
    if len(descriptions) != 0:
        if include_descriptions:
            ticks = np.hstack((x+0.3, x+0.6))
        else:
            ticks = x + 0.4
        ticks.sort()
        labels = []
        for i in this_df.index:
            desc = descriptions.loc[i, 'description']
            if include_descriptions:
                labels.append('('+ desc.replace('\\n', '\n') + ')')

            this_label = i
            this_label = this_label.replace('v1.4_10yrs.db', '')
            this_label = this_label.rstrip('_')
            if plot_type == 'paper':
                this_label = this_label.replace('_', '\_')
            short_name = descriptions.loc[i, 'short']
            if use_short_names:
                this_label = short_name
            labels.append(this_label)
      
    else:
        ticks = x+0.4
        labels = this_df.index
        
    if orientation == 'vertical':
        if plot_type == 'paper' and include_descriptions:
            for i in range(1, len(labels), 2):
                labels[i] = '\\textbf{' + labels[i]+'}'
        # print(labels)
        plt.yticks(ticks=ticks, labels=labels, rotation=0)
        yticks = plt.gca().get_yticklabels()
        for ytck in yticks[1::2]:
            ytck.set_weight('bold')
        plt.gca().tick_params(axis='y', which='both', length=0)
        plt.xlabel(ylabel)
        
        plt.xlim(lims)

        new_ylim = [index_lims[0]*expand_ylim_by[0], 
                    index_lims[1]*expand_ylim_by[1]]
        line_ylim = [-0.1, len(this_df.index)-0.1]
                    
        plt.ylim(new_ylim)

        plt.plot([0,0], line_ylim, 'k', alpha=0.3, zorder=0)

    else:
        plt.xticks(ticks=ticks, labels=labels, rotation=90)
#         xticks = plt.gca().get_xticklabels()
#         for xtck in xticks[1::2]:
#             xtck.set_weight('bold')
        plt.ylabel(ylabel)
        plt.ylim(lims)
        plt.xlim(index_lims)

    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(handles[::-1], labels[::-1], loc=legend_location)

    plt.tight_layout()
    if len(fname) != 0:
        fig.savefig(fname)

def single_bar_graph(df, metrics_display_names, probes=[],
                     plot_type='paper',figsize=(10,15), ylabel='',
                     savefig_path='', save_plot=False,
                     orientation='horizontal',
                     print_title=False, expand_ylim_by=[1.0, 1.0]):

    fig = plt.figure(figsize=figsize)
    wd = 0.8/(len(df.columns))
    x = np.arange(len(df.index))
    x1extra = np.arange(len(df.index)+1)
    k=0
    ticks = []
    labels = []

    

    if len(probes) == 0:
        probes = list(probes_translations.keys())

    for probe in probes:
        c = probes_translations[probe]
        color = colors[probe]
        if c in df.columns:
            lbl = metrics_display_names[c]
            if plot_type == 'paper':
                lbl = lbl.replace('<', '$\\textgreater$')
                lbl = lbl.replace('>', '$\\textless$')
                lbl = lbl.replace('_{', '$_{')
                lbl = lbl.replace('}', '}$')
            # if len(list(baseline_values.keys())) != 0:
            #     lbl += ' (%.3g)' %baseline_values[c]

            if orientation == 'vertical':
                plt.barh(x+k*wd, df[c], height=wd, align='edge', 
                         label=lbl, color=color)
            else:
                plt.bar(x+k*wd, df[c], width=wd, align='edge', 
                         label=lbl, color=color)
            ticks.append((k+0.5)*wd)
            if orientation == 'vertical':
                labels.append(metrics_display_names.loc[c])
            else:
                labels.append(probe.upper())

            k+=1

    index_lims = [0, k*wd]

    if orientation == 'vertical':
        plt.plot([0,0], index_lims, 'k', alpha=0.3)
        plt.xlabel(ylabel)
        plt.yticks(ticks=ticks, labels=labels, rotation=0)
    else:
        plt.plot(index_lims, [0,0],  'k', alpha=0.3)
        plt.ylabel(ylabel)
        plt.xticks(ticks=ticks, labels=labels, rotation=0)
        plt.legend()

    if print_title:
        ttl = df.index[0]
        if plot_type == 'paper':
            ttl = ttl.replace('_', '\_')
            ttl = ttl.replace('.db', '')
        plt.title(ttl)
    plt.tight_layout()

    # new_ylim = [index_lims[0]*expand_ylim_by[0], 
    #                 index_lims[1]*expand_ylim_by[1]]
    ylims = plt.gca().get_ylim()
    new_ylim = [ylims[0]*expand_ylim_by[0], 
                    ylims[1]*expand_ylim_by[1]]
    plt.ylim(new_ylim)

    if len(savefig_path) != 0 and save_plot:
        plt.savefig(savefig_path)


def comparison_plot(metrics_rel, metrics_display_names, misc_metrics, x_axis, 
                    probes=[], themes={}, 
                    ylabel='', xlabel='', savefig_path='', save_plot=False,
                    smoothing='none', nbins=20,
                    plot_baseline=True, expand_ylim_by=[1.1, 1.1],
                    annotations={}):
    """
    probes: list
        Possible values: ['fom', 'lss', 'wl', 'sl', 'sn', 'kn', 'pz']
    """

    sim_annotation_numbers = {
        'footprint_big_sky_dustv1.5_10yrs': 1,
        'footprint_newAv1.5_10yrs' : 2,
        'wfd_depth_scale0.65_noddf_v1.5_10yrs' : 3,
        'wfd_depth_scale0.99_noddf_v1.5_10yrs' : 4,
        'baseline_2snaps_v1.5_10yrs' : 5,
        'footprint_newBv1.5_10yrs' : 6,
        'footprint_bluer_footprintv1.5_10yrs' : 7,
        'bulges_bs_v1.5_10yrs' : 8,
        'bulges_cadence_i_heavy_v1.5_10yrs' : 9,
        'wfd_depth_scale0.99_v1.5_10yrs' : 10,
        'short_exp_5ns_5expt_v1.5_10yrs' : 11,
        'dcr_nham2_ugr_v1.5_10yrs' : 12,
        'dcr_nham1_ugri_v1.5_10yrs' : 13 
    }
    

    themes_translations = {0:t0, 1:t1}
    
    if len(list(themes.keys())) == 0:
        for p in probes:
            themes[p] = 0

    for i in misc_metrics.index:
        if np.isnan(misc_metrics.loc[i, x_axis]):
            misc_metrics = misc_metrics.drop(i)
            metrics_rel = metrics_rel.drop(i)

    xvals = misc_metrics.loc[metrics_rel.index, x_axis].astype('float')
    inds = np.argsort(xvals)

    plt.figure()
    

    ylim = [0,0]

    for i in range(len(probes)):
        metric = probes_translations[probes[i]]

        yvals = metrics_rel.loc[metrics_rel.index, metric].astype('float')

        x = xvals[inds].values
        y = yvals[inds].values

        unique_x = np.unique(x)
        unique_y = np.zeros(len(unique_x))

        for k in range(len(unique_x)):
            msk = x == unique_x[k]
            if sum(msk) == 1:
                unique_y[k] = y[msk]
            else:
                unique_y[k] = y[msk].mean()

        x = unique_x
        y = unique_y

        if smoothing == 'none':
            new_x = x
            new_y = y
            yerr = np.zeros(len(new_x))

        elif smoothing == 'interpolate':
            # Ok this just doesn't work, data is too noisy
            

            f = interp1d(unique_x, unique_y, kind='cubic')
            new_x = np.linspace(x.min(), x.max(), 100)
            new_y = f(new_x)

        elif smoothing == 'gp':
            kernel = george.kernels.ExpSquaredKernel((x.max()-x.min())/100)
            gp = george.GP(kernel, fit_mean=True, fit_white_noise=True)
            gp.compute(x, y.std())

            def neg_ln_like(p):
                gp.set_parameter_vector(p)
                return -gp.log_likelihood(y)

            def grad_neg_ln_like(p):
                gp.set_parameter_vector(p)
                return -gp.grad_log_likelihood(y)

            result = minimize(neg_ln_like, gp.get_parameter_vector(),
                                      jac=grad_neg_ln_like, method="L-BFGS-B",
                                      tol = 1e-5)
            gp.set_parameter_vector(result.x)

            pred_x = np.linspace(x.min(), x.max(), 100)
            pred_y, y_var = gp.predict(y, pred_x, return_var=True)
            new_x = pred_x
            new_y = pred_y

        elif smoothing == 'binning':
            xgrid = np.linspace(x.min(), x.max(), nbins)
            x_delt = xgrid[1] - xgrid[0]
            new_x = []
            new_y = []
            yerr = []

            for k in range(len(xgrid)):
                msk = (x>=xgrid[k]) & (x<xgrid[k] + x_delt)
                if msk.sum() == 0:
                    pass
                elif msk.sum() == 1:
                    new_x.append(x[msk][0])
                    new_y.append(y[msk][0])
                    yerr.append(0)
                else:
                    new_x.append(x[msk].mean())
                    new_y.append(y[msk].mean())
                    yerr.append(y[msk].std())

            new_x = np.array(new_x)
            new_y = np.array(new_y)
            yerr = np.array(yerr)
            

        theme = themes_translations[themes[probes[i]]]
        plt.errorbar(new_x, new_y, yerr=yerr, marker=None,
                  lw=theme['lw'], alpha=theme['alpha'], 
                 color=colors[probes[i]])
        plt.plot(new_x, new_y, marker=markers[probes[i]],
                 label=metrics_display_names.loc[metric], lw=0,
                 linestyle=None, alpha=theme['alpha'], color=colors[probes[i]])

        if (new_y+yerr).max() > ylim[1]:
            ylim[1] = (new_y+yerr).max()
        if (new_y-yerr).min() < ylim[0]:
            ylim[0] = (new_y-yerr).min()
        

    if len(annotations) != 0:
        anno_x = []
        anno_y = []
        anno_sim = []
        anno_num = []
        for k in annotations.keys():
            which_probe = annotations[k][0]
            metric = probes_translations[which_probe]
            yshift = -0.06 + annotations[k][1]
            xtext = misc_metrics.loc[k, x_axis].astype('float')
            ytext = metrics_rel.loc[k, metric].astype('float') + yshift
            anno_x.append(xtext)
            anno_y.append(ytext)
            anno_sim.append(k)
            anno_num.append(sim_annotation_numbers[k])

        sorted_inds = np.argsort(anno_x)
        counter = 1
        print_str = 'Annotations: '
        for ind in sorted_inds:
            plt.text(anno_x[ind], anno_y[ind], (str)(counter), color='k')
            # plt.text(anno_x[ind], anno_y[ind], (str)(anno_num[ind]), color='k')
            sim_name = anno_sim[ind]
            sim_name = sim_name.replace('_', '\_')
            print_str += '%d-\\texttt{%s}, ' %(counter, sim_name)
            # print_str += '%d-\\texttt{%s}, ' %(anno_num[ind], sim_name)
            counter += 1
        print(print_str[:-2]) 

    # handles, labels = plt.gca().get_legend_handles_labels()
    # # remove the errorbars
    # print(handles[0][0])
    # handles = [h[0] for h in handles]
    # print(handles[0])
 
    # plt.legend(handles, labels, numpoints=1, loc='best')
    plt.legend(loc='best')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    

    ax = plt.gca()
    new_ylim = [ylim[0]*expand_ylim_by[0], ylim[1]*expand_ylim_by[1]]
    new_xlim = ax.get_xlim()
    
   
    zo = 0

    plt.plot(new_xlim, [0, 0], '--k', alpha=0.5, zorder=zo)

    if plot_baseline:

        base_val = float(misc_metrics.loc[read_metrics.baseline, x_axis])
        plt.plot([base_val, base_val], new_ylim, '--k', lw=1, alpha=0.3,
                 zorder=zo)
        plt.text(base_val, new_ylim[1]*0.8, 'baseline', 
                alpha=0.3, zorder=zo)


    plt.xlim(new_xlim)
    plt.ylim(new_ylim)
    plt.tight_layout()
    
    if len(savefig_path) != 0 and save_plot:
        plt.savefig(savefig_path)
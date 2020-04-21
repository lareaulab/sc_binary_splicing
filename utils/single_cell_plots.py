import numpy as np
import pandas as pd
import subprocess as sp
import os
from matplotlib import pyplot as plt
from scipy.stats import pearsonr
from matplotlib.gridspec import GridSpec
import matplotlib.cm as cm
import seaborn as sns
from IPython.core.pylabtools import figsize
from pylab import *

import matplotlib.colors as mcolors


import matplotlib as mpl
import numpy as np
from matplotlib import pyplot as plt

mpl.rcParams["mathtext.fontset"] = "stix"


def plot_pca(pca_tab, psi_tab, mrna_tab, reads_tab, cj, event, mrna_min, plot_dir, cmap='viridis', 
             filter_cells=True, print_pseudotime=False, vmin=0, vmax=1, alpha_missing=0.2, width = 8.5):
    
#     figsize(10,8)
    
    #pass_cells = psi_tab.columns[(mrna_tab.loc[event] > mrna_min) & (rpm_tab.loc[event] >= rpm_tab.mean(axis=1).quantile(0.1))]
    
    
    if filter_cells:
        read_min = mrna_min*cj*(1+psi_tab.loc[event])

        pass_cells = psi_tab.columns[(mrna_tab.loc[event] > mrna_min) & (reads_tab.loc[event] >= read_min)]
        
    else:
        pass_cells = psi_tab.columns[(mrna_tab.loc[event] > -1) & (reads_tab.loc[event] >= -1)]


    figsize(width, 5)
    fig = plt.figure()
    ax  = plt.subplot(1,1,1)
    ax.scatter(pca_tab.PC1, pca_tab.PC2, c='gray', s=75, edgecolors='none', alpha=alpha_missing)
    sc = ax.scatter(pca_tab.loc[pass_cells, 'PC1'], pca_tab.loc[pass_cells, 'PC2'], 
                c=psi_tab.loc[event, pass_cells], s=100, edgecolors='none', vmin=vmin, vmax=vmax, cmap=cmap)
    
    
    
    

    ax.tick_params(labelsize=28, length=5)
    cb = plt.colorbar(sc, shrink = 0.8)
    cb.set_label(label='$\Psi$',size=28)
    cb.ax.tick_params(labelsize=28, length=5)
    cb.outline.set_visible(False)
    

    #plt.plot(chen_pca.line_1, -chen_pca.line_2, c='red',linewidth=3, label='lineage')
    
    ax.set_xlim([-65, 102])
    ax.set_ylim([-55, 63])
    
    
#     cb = plt.colorbar(sc, shrink = 0.8)
#     cb.set_label(label='$\Psi$',size=28)
#     cb.ax.tick_params(labelsize=28, length=5)
#     cb.outline.set_visible(False)
    
    #plt.legend(frameon=False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    
    
    ax.tick_params(labelsize=28, length=5)

    plt.ylabel('PC2', fontsize=28)
    plt.xlabel('PC1', fontsize=28)
    plt.title(event.split('_')[0], fontsize=28)
    
    if filter_cells:

        plt.savefig(plot_dir+event+'.pca.svg', bbox_inches='tight', transparent=True)
        plt.savefig(plot_dir+event+'.pca.pdf', bbox_inches='tight', transparent=True)
        plt.savefig(plot_dir+event+'.pca.png', dpi=300, bbox_inches='tight', transparent=True)
        
    else:
        
        plt.savefig(plot_dir+event+'.pca.no_filter.svg', bbox_inches='tight', transparent=True)
        plt.savefig(plot_dir+event+'.pca.no_filter.pdf', bbox_inches='tight', transparent=True)
        plt.savefig(plot_dir+event+'.pca.no_filter.png', dpi=300, bbox_inches='tight', transparent=True)

        
    
    
    #plt.show()
    
    if print_pseudotime:
    
        ####figsize(6, 4)
        figsize(6, 4)
        fig = plt.figure()
        ax  = plt.subplot(1,1,1)
        ax.scatter(pca_tab.pseudotime, psi_tab.loc[event], c='skyblue', s=75, edgecolors='none', alpha=0.5)
        ax.scatter(pca_tab.loc[pass_cells, 'pseudotime'], psi_tab.loc[event, pass_cells], 
                    c='navy', s=150, edgecolors='none')

        ax.set_ylim((-0.05,1.05))

        #plt.plot(chen_pca.line_1, -chen_pca.line_2, c='red',linewidth=3, label='lineage')
        #plt.legend(frameon=False)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)

        ax.tick_params(labelsize=28, length=5)

        plt.ylabel('$\Psi$', fontsize=28)
        plt.xlabel('pseudotime', fontsize=28)
        plt.title(event.split('_')[0], fontsize=28)

        if filter_cells:

            plt.savefig(plot_dir+event+'.pseudotime.svg', bbox_inches='tight', transparent=True)
            plt.savefig(plot_dir+event+'.pseudotime.pdf', bbox_inches='tight', transparent=True)
            plt.savefig(plot_dir+event+'.pseudotime.png', dpi=300, bbox_inches='tight', transparent=True)

        else:
            plt.savefig(plot_dir+event+'.pseudotime.no_filter.svg', bbox_inches='tight', transparent=True)
            plt.savefig(plot_dir+event+'.pseudotime.no_filter.pdf', bbox_inches='tight', transparent=True)
            plt.savefig(plot_dir+event+'.pseudotime.no_filter.png', dpi=300, bbox_inches='tight', transparent=True)
    
    plt.show()
    

def get_bins_table(PSI_table, total_counts, nbins = 25, int_min = 0.2):
    
    '''
    This function gets a sorted heatmap of PSI distribution, ready to plot.
    output
    '''
    
    steps = 1/nbins

    expression_sort = total_counts.mean(axis=1).sort_values().index

    bins_data = pd.DataFrame()
    #psi_bins = np.arange(0, 1 + 2*(steps), steps)
    psi_bins = np.append(np.arange(0, 1, steps), 1) ####
    # np.append(np.arange(0, 1, 1/11), 1)
    for i in range(len(psi_bins)-1):
        lower = psi_bins[i]
        upper = psi_bins[i+1]
        #bin_counts = ((PSI_table >= lower) & (PSI_table < upper)).sum(axis = 1)
        bin_counts = ((PSI_table >= lower) & (PSI_table <= upper)).sum(axis = 1) ####
        max_cells = (len(PSI_table.columns)-PSI_table.isna().sum(axis=1))
        
        #### This box
        if upper == 1:
            bin_name = 'bin_1.0'
        else:
            bin_name = 'bin_' + str(round(lower, 2))
        ####
            
        #bins_data['bin_' + str(round(lower, 2))] = bin_counts/max_cells
        bins_data[bin_name] = bin_counts/max_cells ####
        

    bin_norm = bins_data.div(bins_data.max(axis=1), axis=0)
    hist_complete = bin_norm.loc[expression_sort]
    #hist_complete.columns = [str(round(x, 2)) for x in np.arange(0, 1 + steps, steps)]
    col_names = [x.split('_')[1] for x in bins_data.columns] ####
    if nbins%2 == 1:
        col_names[nbins//2] = '0.5' ####
    #hist_complete.columns = [x.split('_')[1] for x in bins_data.columns] ####
    hist_complete.columns = [col_names] ####
    

    filtro = (PSI_table.mean(axis=1) >=int_min) & (PSI_table.mean(axis=1) <=(1-int_min))
    expression_fil = total_counts.loc[filtro]
    PSI_fil = PSI_table.loc[filtro]


    expression_sort = expression_fil.mean(axis=1).sort_values().index

    bins_data = pd.DataFrame()
    #psi_bins = np.arange(0, 1 + (2*steps), steps)
    psi_bins = np.append(np.arange(0, 1, steps), 1) ####
    for i in range(len(psi_bins)-1):
        lower = psi_bins[i]
        upper = psi_bins[i+1]
        #bin_counts = ((PSI_fil >= lower) & (PSI_fil < upper)).sum(axis = 1)
        bin_counts = ((PSI_fil >= lower) & (PSI_fil <= upper)).sum(axis = 1)
        max_cells = (len(PSI_fil.columns)-PSI_fil.isna().sum(axis=1))
        if upper == 1:
            bin_name = 'bin_1.0'
        else:
            bin_name = 'bin_' + str(round(lower, 2))
            
        #bins_data['bin_' + str(round(lower, 2))] = bin_counts/max_cells
        bins_data[bin_name] = bin_counts/max_cells
        
    bin_norm = bins_data.div(bins_data.max(axis=1), axis=0)
    hist_intermediate = bin_norm.loc[expression_sort]
    #hist_intermediate.columns = [str(round(x, 2)) for x in np.arange(0, 1+steps, steps)]
    col_names = [x.split('_')[1] for x in bins_data.columns] ####
    if nbins%2 == 1:
        col_names[nbins//2] = '0.5' ####
    #hist_intermediate.columns = [x.split('_')[1] for x in bins_data.columns] ####
    hist_intermediate.columns = [col_names] ####
    
    return hist_complete, hist_intermediate


def get_bins_table2(PSI_table, total_counts, nbins = 25, int_min = 0.2):
    
    '''
    This function gets a sorted heatmap of PSI distribution, ready to plot.
    output
    '''
    
    steps = 1/nbins
    
    expression_sort = PSI_table.mean(axis=1).sort_values().index

    bins_data = pd.DataFrame()
    #psi_bins = np.arange(0, 1 + 2*(steps), steps)
    psi_bins = np.append(np.arange(0, 1, steps), 1)
    for i in range(len(psi_bins)-1):
        lower = psi_bins[i]
        upper = psi_bins[i+1]
        #bin_counts = ((PSI_table >= lower) & (PSI_table < upper)).sum(axis = 1)
        bin_counts = ((PSI_table >= lower) & (PSI_table <= upper)).sum(axis = 1) ####
        max_cells = (len(PSI_table.columns)-PSI_table.isna().sum(axis=1))
        
        #### This box
        if upper == 1:
            bin_name = 'bin_1.0'
        else:
            bin_name = 'bin_' + str(round(lower, 2))
        ####
            
        #bins_data['bin_' + str(round(lower, 2))] = bin_counts/max_cells
        bins_data[bin_name] = bin_counts/max_cells ####
        

    bin_norm = bins_data.div(bins_data.max(axis=1), axis=0)
    hist_complete = bin_norm.loc[expression_sort]
    #hist_complete.columns = [str(round(x, 2)) for x in np.arange(0, 1 + steps, steps)]
    col_names = [x.split('_')[1] for x in bins_data.columns] ####
    if nbins%2 == 1:
        col_names[nbins//2] = '0.5' ####
    #hist_complete.columns = [x.split('_')[1] for x in bins_data.columns] ####
    hist_complete.columns = [col_names] ####
    

    filtro = (PSI_table.mean(axis=1) >=int_min) & (PSI_table.mean(axis=1) <=(1-int_min))
    expression_fil = total_counts.loc[filtro]
    PSI_fil = PSI_table.loc[filtro]


    expression_sort = PSI_fil.mean(axis=1).sort_values().index

    bins_data = pd.DataFrame()
    #psi_bins = np.arange(0, 1 + (2*steps), steps)
    psi_bins = np.append(np.arange(0, 1, steps), 1)
    for i in range(len(psi_bins)-1):
        lower = psi_bins[i]
        upper = psi_bins[i+1]
        #bin_counts = ((PSI_fil >= lower) & (PSI_fil < upper)).sum(axis = 1)
        bin_counts = ((PSI_fil >= lower) & (PSI_fil <= upper)).sum(axis = 1)
        max_cells = (len(PSI_fil.columns)-PSI_fil.isna().sum(axis=1))
        if upper == 1:
            bin_name = 'bin_1.0'
        else:
            bin_name = 'bin_' + str(round(lower, 2))
            
        #bins_data['bin_' + str(round(lower, 2))] = bin_counts/max_cells
        bins_data[bin_name] = bin_counts/max_cells
        
    bin_norm = bins_data.div(bins_data.max(axis=1), axis=0)
    hist_intermediate = bin_norm.loc[expression_sort]
    #hist_intermediate.columns = [str(round(x, 2)) for x in np.arange(0, 1+steps, steps)]
    col_names = [x.split('_')[1] for x in bins_data.columns] ####
    if nbins%2 == 1:
        col_names[nbins//2] = '0.5' ####
    #hist_intermediate.columns = [x.split('_')[1] for x in bins_data.columns] ####
    hist_intermediate.columns = [col_names] ####
    
    return hist_complete, hist_intermediate


def plot_cell_expression_v_binary(PSI_table, total_counts, dset_name, save_name, r_pos = False, 
                                  xlocations=False, xtags=False, xlabel = 'total SJ reads in cell (log10)', 
                                  xlim = False, plot_dir = 'plots/', just_show=False, binary_def = 0):
    
    figsize(7, 6)
    
    if not os.path.isdir(plot_dir):
        os.makedirs(plot_dir)
            
    observed = np.sum(total_counts).index[np.sum(total_counts) > 0]
    
    total_counts = total_counts[observed]
    PSI_table = PSI_table[observed]
    
    counts_mean = np.sum(total_counts)
    
    psi_mean = (np.sum((PSI_table >= (1-binary_def)) | (PSI_table <= (0+binary_def))))/np.sum(~PSI_table.isnull())
     
    good_psi = psi_mean.dropna()
        
    pr = pearsonr(np.log10(np.array(counts_mean.loc[good_psi.index])), good_psi)[0]
    fig = plt.figure()
    
    #fig.suptitle(dset_name, fontsize=28, y=0.95, x=0.43)

    gs = GridSpec(4,4)
    
    gs.update(wspace=0.05, hspace=0.05)

    ax_joint = fig.add_subplot(gs[1:4,0:3])
    
    ax_marg_x = fig.add_subplot(gs[0,0:3])
    ax_marg_x.spines["top"].set_visible(False)
    ax_marg_x.spines["bottom"].set_visible(False)
    ax_marg_x.spines["left"].set_visible(False)
    ax_marg_x.spines["right"].set_visible(False)
    ax_marg_x.xaxis.set_ticks_position('none') 
    ax_marg_x.yaxis.set_ticks_position('none') 
    
    ax_marg_y = fig.add_subplot(gs[1:4,3])
    ax_marg_y.spines["top"].set_visible(False)
    ax_marg_y.spines["bottom"].set_visible(False)
    ax_marg_y.spines["left"].set_visible(False)
    ax_marg_y.spines["right"].set_visible(False)
    ax_marg_y.xaxis.set_ticks_position('none') 
    ax_marg_y.yaxis.set_ticks_position('none') 
    
    
    
    ax_joint.scatter(np.log10(np.array(counts_mean)), psi_mean, c='navy', alpha=0.5, s=50, edgecolors='none')
    ax_joint.tick_params(labelsize=28, length=5)
    ax_marg_x.hist(np.log10(np.array(counts_mean)), color='skyblue')
    ax_marg_y.hist(psi_mean.dropna(),orientation="horizontal", color='skyblue')
    ax_marg_y.set_ylim((-0.05,1.05))
    ax_joint.set_ylim((-0.05,1.05))
    
    if xlim:
        ax_marg_x.set_xlim(xlim)
        ax_joint.set_xlim(xlim)
        
#     else:
#         ax_marg_x.set_xlim((np.min(np.log10(np.array(counts_mean))), np.max(np.log10(np.array(counts_mean)))))
#         ax_joint.set_xlim((np.min(np.log10(np.array(counts_mean))), np.max(np.log10(np.array(counts_mean)))))

    # Turn off tick labels on marginals
    plt.setp(ax_marg_x.get_xticklabels(), visible=False)
    plt.setp(ax_marg_x.get_yticklabels(), visible=False)
    plt.setp(ax_marg_y.get_xticklabels(), visible=False)
    plt.setp(ax_marg_y.get_yticklabels(), visible=False)
    
    if not r_pos:
        rpos_max = np.log10(np.array(counts_mean)).max()
        rpos_min = np.log10(np.array(counts_mean)).min()
        r_pos = rpos_min + 0.1*(rpos_max - rpos_min)
    
    ax_joint.text(r_pos, 0.2, 'R = '+str(round(pr,3)), fontsize=28, verticalalignment='top')

    # Set labels on joint
    ax_joint.set_xlabel(xlabel, fontsize=28)
    ax_joint.set_ylabel(r'% exons with binary $\Psi$', fontsize=28)
    
    ax_joint.spines["right"].set_visible(False)
    ax_joint.spines["top"].set_visible(False)
#     ax_joint.spines['left'].set_linewidth(2)
#     ax_joint.spines['bottom'].set_linewidth(2)

#     ax_joint.xaxis.set_tick_params(width=2)
#     ax_joint.yaxis.set_tick_params(width=2)



    
    
    if xlocations:
        
        ax_joint.set_xticks([np.log10(x) for x in xlocations])
        ax_joint.set_xticklabels(xtags)
    
#         xlocs = np.array([np.log10(x) for x in xlocations])
#         if xtags:
#             xnewLabels = xtags
#         else:
#             xnewLabels = np.array([str(x) for x in xlocations])   
        
#         plt.sca(ax_joint)

#         plt.xticks(xlocs, xnewLabels)

    
    if not just_show:
        plt.savefig(plot_dir + '/' + save_name + '_cell_expression_v_binary.svg', bbox_inches='tight', transparent=True)
        plt.savefig(plot_dir + '/' + save_name + '_cell_expression_v_binary.pdf', bbox_inches='tight', transparent=True)
        plt.savefig(plot_dir + '/' + save_name + '_cell_expression_v_binary.png', bbox_inches='tight', transparent=True)
        
    else:
        plt.show()
        
    return pr
        
        
        

def plot_ase_expression_v_binary(PSI_table, total_counts, dset_name, save_name, 
                                 xlocations=False, xtags=False, xlabel = 'total SJ reads in cell (log10 + 1)', 
                                 xlim = False, plot_dir = 'plots/', just_show=False, binary_def = 0):
    figsize(7, 6)
    if not os.path.isdir(plot_dir):
        os.makedirs(plot_dir)
        
        
    observed = np.sum(total_counts).index[np.sum(total_counts) > 0]
    
    total_counts = total_counts[observed]
    PSI_table = PSI_table[observed]
    
    counts_mean = []
    psi_mean = []
    
    y = ( (PSI_table<=(0+binary_def)).sum(axis=1) + (PSI_table>=(1-binary_def)).sum(axis=1) )/(~PSI_table.isnull()).sum(axis=1)
    
    fig = plt.figure()
    
    #fig.suptitle(dset_name, fontsize=28, y=0.95, x=0.43)

    gs = GridSpec(4,4)
    
    gs.update(wspace=0.05, hspace=0.05)

    ax_joint = fig.add_subplot(gs[1:4,0:3])
    ax_marg_x = fig.add_subplot(gs[0,0:3])
    ax_marg_x.spines["top"].set_visible(False)
    ax_marg_x.spines["bottom"].set_visible(False)
    ax_marg_x.spines["left"].set_visible(False)
    ax_marg_x.spines["right"].set_visible(False)
    ax_marg_x.xaxis.set_ticks_position('none') 
    ax_marg_x.yaxis.set_ticks_position('none') 
    
    ax_marg_y = fig.add_subplot(gs[1:4,3])
    ax_marg_y.spines["top"].set_visible(False)
    ax_marg_y.spines["bottom"].set_visible(False)
    ax_marg_y.spines["left"].set_visible(False)
    ax_marg_y.spines["right"].set_visible(False)
    ax_marg_y.xaxis.set_ticks_position('none') 
    ax_marg_y.yaxis.set_ticks_position('none') 
    
    pr = pearsonr(np.log10(total_counts.mean(axis=1)+1), y)[0]
    
    
    colors1 = plt.cm.viridis_r(np.linspace(0., 1, 128))
    colors2 = plt.cm.viridis(np.linspace(0, 1, 128))

    # combine them and build a new colormap
    colors = np.vstack((colors1, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

#     plt.pcolor(data, cmap=mymap)
#     plt.colorbar()
#     plt.show()
    
    
    sc = ax_joint.scatter(np.log10(total_counts.mean(axis=1)+1), y, c=PSI_table.mean(axis=1), s=50, edgecolors='none',
                         cmap=mymap, vmin = 0, vmax=1)
    ax_joint.tick_params(labelsize=28, length=5)
    cb = plt.colorbar(sc, shrink = 1.8)
    cb.set_label(label='average exon $\Psi$',size=28)
    cb.ax.tick_params(labelsize=28, length=5)
    cb.outline.set_visible(False)
    
    ax_marg_x.hist(np.log10(total_counts.mean(axis=1)+1), color='navy')
    ax_marg_y.hist(y.dropna(),orientation="horizontal", color='gold')
    ax_marg_y.set_ylim((-0.05,1.05))
    ax_joint.set_ylim((-0.05,1.05))

    # Turn off tick labels on marginals
    plt.setp(ax_marg_x.get_xticklabels(), visible=False)
    plt.setp(ax_marg_x.get_yticklabels(), visible=False)
    plt.setp(ax_marg_y.get_xticklabels(), visible=False)
    plt.setp(ax_marg_y.get_yticklabels(), visible=False)
    
    # Set labels on joint
    ax_joint.set_xlabel(xlabel, fontsize=28)
    ax_joint.set_ylabel(r'% cells with binary $\Psi$', fontsize=28)
    
    ax_joint.spines["right"].set_visible(False)
    ax_joint.spines["top"].set_visible(False)
#     ax_joint.spines['left'].set_linewidth(2)
#     ax_joint.spines['bottom'].set_linewidth(2)

#     ax_joint.xaxis.set_tick_params(width=2)
#     ax_joint.yaxis.set_tick_params(width=2)
    
    if xlocations:
    
        xlocs = np.array([np.log10(x+1) for x in xlocations])
        if xtags:
            xnewLabels = xtags
        else:
            xnewLabels = np.array([str(x) for x in xlocations])   
        
        plt.sca(ax_joint)

        plt.xticks(xlocs, xnewLabels)


    plt.savefig(plot_dir + '/' + save_name + '_ase_expression_v_binary.svg', bbox_inches='tight', transparent=True)
    plt.savefig(plot_dir + '/' + save_name + '_ase_expression_v_binary.pdf', bbox_inches='tight', transparent=True)
    plt.savefig(plot_dir + '/' + save_name + '_ase_expression_v_binary.png', bbox_inches='tight', transparent=True)
    
    return pr
    
    
def plot_all_expression_v_psi(PSI_table, total_counts, dset_name, save_name, alpha=0.02, plot_dir = 'plots/', just_show=False):
    
    if not os.path.isdir(plot_dir):
        os.makedirs(plot_dir)
    
    #counts_mean = []
    #psi_mean = []
    
    #y = ( (PSI_table==0).sum(axis=1) + (PSI_table==1).sum(axis=1) )/(~PSI_table.isnull()).sum(axis=1)
    
    counts = np.log1p(np.array(total_counts).flatten())
    
    psi = np.abs(0.5 - np.array(PSI_table).flatten())
    
    fig = plt.figure()
    
    fig.suptitle(dset_name, fontsize=28, y=0.95, x=0.43)

    gs = GridSpec(4,4)

    ax_joint = fig.add_subplot(gs[1:4,0:3])
    ax_marg_x = fig.add_subplot(gs[0,0:3])
    ax_marg_y = fig.add_subplot(gs[1:4,3])
    

    ax_joint.scatter(counts, psi, s=100, edgecolors='none', alpha=alpha, c='navy')
    #plt.colorbar(sc, label='abs(0.5-$\mu(\Psi)$)')
    
    ax_marg_x.hist(counts, color='navy')
    ax_marg_y.hist(psi,orientation="horizontal", color='gold')
    #ax_marg_y.set_ylim((0,1))
    #ax_joint.set_ylim((0,1))

    # Turn off tick labels on marginals
    plt.setp(ax_marg_x.get_xticklabels(), visible=False)
    plt.setp(ax_marg_y.get_yticklabels(), visible=False)
    
    # Set labels on joint
    ax_joint.set_xlabel(r'$\mu$' +' (log1p(counts)) per ASE per cell', fontsize=18)
    ax_joint.set_ylabel(r'abs(0.5-$\mu(\Psi)$)', fontsize=18)

    plt.savefig(plot_dir + '/' + save_name + '_all_expression_v_psi.svg', bbox_inches='tight', transparent=True)
    plt.savefig(plot_dir + '/' + save_name + '_all_expression_v_psi.pdf', bbox_inches='tight', transparent=True)
    plt.savefig(plot_dir + '/' + save_name + '_all_expression_v_psi.png', bbox_inches='tight', transparent=True)
    
    
def plot_ase_expression_v_missing(PSI_table, total_counts, dset_name, save_name, plot_dir='plots/', just_show=False):
    
    
    
    if not os.path.isdir(plot_dir):
        os.makedirs(plot_dir)
        
    counts_mean = []
    psi_mean = []
    
    y = PSI_table.isnull().mean(axis=1)

    fig = plt.figure()
    
    fig.suptitle(dset_name, fontsize=28, y=0.95, x=0.43)

    gs = GridSpec(4,4)

    ax_joint = fig.add_subplot(gs[1:4,0:3])
    ax_marg_x = fig.add_subplot(gs[0,0:3])
    ax_marg_y = fig.add_subplot(gs[1:4,3])
    

    ax_joint.scatter(np.log1p(total_counts.mean(axis=1)), y, c='navy', alpha=0.5, s=100, edgecolors='none')
    
    ax_marg_x.hist(np.log1p(total_counts.mean(axis=1)), color='skyblue')
    ax_marg_y.hist(y,orientation="horizontal", color='skyblue')
    ax_marg_y.set_ylim((0,1))
    ax_joint.set_ylim((0,1))

    # Turn off tick labels on marginals
    plt.setp(ax_marg_x.get_xticklabels(), visible=False)
    plt.setp(ax_marg_y.get_yticklabels(), visible=False)
    
    #ax_joint.text(x, 0.2, 'R = '+str(round(pr,3)), fontsize=20, verticalalignment='top')

    # Set labels on joint
    ax_joint.set_xlabel(r'$\mu$' +' (log1p(counts)) per ASE', fontsize=18)
    ax_joint.set_ylabel(r'Proportion of cells missing ASE', fontsize=18)

    if not just_show:
        plt.savefig(plot_dir + '/' + save_name + '_ase_expression_v_missing.svg', bbox_inches='tight', transparent=True)
        plt.savefig(plot_dir + '/' + save_name + '_ase_expression_v_missing.pdf', bbox_inches='tight', transparent=True)
        plt.savefig(plot_dir + '/' + save_name + '_ase_expression_v_missing.png', bbox_inches='tight', transparent=True)
        
    else:
        plt.show()
        
    
def plot_cell_expression_v_missing(PSI_table, total_counts, dset_name, save_name, plot_dir = 'plots/', just_show=False):
    
    if not os.path.isdir(plot_dir):
        os.makedirs(plot_dir)
        
    counts_mean = []
    psi_mean = []
    
    y = PSI_table.isnull().mean(axis=0)
   
    fig = plt.figure()
    
    fig.suptitle(dset_name, fontsize=28, y=0.95, x=0.43)

    gs = GridSpec(4,4)

    ax_joint = fig.add_subplot(gs[1:4,0:3])
    ax_marg_x = fig.add_subplot(gs[0,0:3])
    ax_marg_y = fig.add_subplot(gs[1:4,3])
    

    ax_joint.scatter(np.log1p(total_counts.mean(axis=0)), y, c='navy', alpha=0.5, s=100, edgecolors='none')
    
    ax_marg_x.hist(np.log1p(total_counts.mean(axis=0)), color='skyblue')
    ax_marg_y.hist(y,orientation="horizontal", color='skyblue')
    ax_marg_y.set_ylim((0,1))
    ax_joint.set_ylim((0,1))

    # Turn off tick labels on marginals
    plt.setp(ax_marg_x.get_xticklabels(), visible=False)
    plt.setp(ax_marg_y.get_yticklabels(), visible=False)
    
    # Set labels on joint
    ax_joint.set_xlabel(r'$\mu$' +' (log1p(counts)) per ASE', fontsize=18)
    ax_joint.set_ylabel(r'Proportion of ASE missing in cell', fontsize=18)

    if not just_show:
        plt.savefig(plot_dir + '/' + save_name + '_cell_expression_v_missing.svg', bbox_inches='tight', transparent=True)
        plt.savefig(plot_dir + '/' + save_name + '_cell_expression_v_missing.pdf', bbox_inches='tight', transparent=True)
        plt.savefig(plot_dir + '/' + save_name + '_cell_expression_v_missing.png', dpi=300, bbox_inches='tight', transparent=True)
        
    else:
        plt.show()
    
    
def plot_compare_dataset_missing(PSI_table1, PSI_table2, dset1, dset2, save_name, plot_dir='plots/', just_show=False):
    
    if not os.path.isdir(plot_dir):
        os.makedirs(plot_dir)
    
    shared = [x for x in PSI_table1.index if x in PSI_table2.index]
    x = PSI_table1.loc[shared].isnull().mean(axis=1)
    y = PSI_table2.loc[shared].isnull().mean(axis=1)
    
    print(len(shared))
    
    fig = plt.figure()
    
    fig.suptitle('ASE missing, ' + dset1 + ' vs ' + dset2, fontsize=28, y=0.95, x=0.43)

    gs = GridSpec(4,4)

    ax_joint = fig.add_subplot(gs[1:4,0:3])
    ax_marg_x = fig.add_subplot(gs[0,0:3])
    ax_marg_y = fig.add_subplot(gs[1:4,3])
    

    ax_joint.scatter(x, y, c='navy', alpha=0.5, s=100, edgecolors='none')
    
    ax_marg_x.hist(x, color='skyblue')
    ax_marg_y.hist(y,orientation="horizontal", color='skyblue')
    ax_marg_y.set_ylim((-0.1,1.1))
    ax_joint.set_ylim((-0.1,1.1))
    
    ax_marg_x.set_xlim((-0.1,1.1))
    ax_joint.set_xlim((-0.1,1.1))

    # Turn off tick labels on marginals
    plt.setp(ax_marg_x.get_xticklabels(), visible=False)
    plt.setp(ax_marg_y.get_yticklabels(), visible=False)
    
    ax_joint.text(0.7, 0.2, 'n = '+str(len(shared)), fontsize=20, verticalalignment='top')

    # Set labels on joint
    ax_joint.set_xlabel('% cells missing ASE in ' + dset1, fontsize=18)
    ax_joint.set_ylabel('% cells missing ASE in ' + dset2, fontsize=18)

    if not just_show:
        plt.savefig(plot_dir + '/' + save_name + '_compare_missing.svg', bbox_inches='tight', transparent=True)
        plt.savefig(plot_dir + '/' + save_name + '_compare_missing.pdf', bbox_inches='tight', transparent=True)
        plt.savefig(plot_dir + '/' + save_name + '_compare_missing.png', bbox_inches='tight', transparent=True)
        
    else:
        plt.show()
        
    
        
def plot_event(PSI_table, total_table, dset, event, xtags, plot_dir = 'plots/', just_show=False):
    
    if not os.path.isdir(plot_dir):
        os.makedirs(plot_dir)
        
    figsize(7, 6)

        
        
    observed = total_table.columns[(total_table.loc[event] > 0)]
    #observed = np.sum(total_table).index[np.sum(total_table) > 0]
    
    total_table = total_table[observed]
    PSI_table = PSI_table[observed]
    
    x = np.log10(total_table.loc[event])
    y = PSI_table.loc[event]
    
    
    fig = plt.figure()
    
    fig.suptitle(event.split('_')[0], fontsize=28, y=0.85, x=0.75)

    gs = GridSpec(4,4)
    
    gs.update(wspace=0.05, hspace=0.05)

    ax_joint = fig.add_subplot(gs[1:4,0:3])
    ax_marg_x = fig.add_subplot(gs[0,0:3])
    ax_marg_y = fig.add_subplot(gs[1:4,3])
    

    ax_joint.scatter(x, y, c='navy', alpha=0.5, s=75, edgecolors='none')
    
    ax_marg_x.hist(x.dropna(), color='skyblue')
    ax_marg_y.hist(y.dropna(),orientation="horizontal", color='skyblue')
    

    # Turn off tick labels on marginals
    plt.setp(ax_marg_x.get_xticklabels(), visible=False)
    plt.setp(ax_marg_y.get_yticklabels(), visible=False)
    
    # Set labels on joint
    ax_joint.set_xlabel('SJ counts', fontsize=28)
    ax_joint.set_ylabel('$\Psi$', fontsize=28)
    
    
    
    
    
    
    
    
    ######
    
    ax_marg_x.spines["top"].set_visible(False)
    ax_marg_x.spines["bottom"].set_visible(False)
    ax_marg_x.spines["left"].set_visible(False)
    ax_marg_x.spines["right"].set_visible(False)
    ax_marg_x.xaxis.set_ticks_position('none') 
    ax_marg_x.yaxis.set_ticks_position('none') 
    
    ax_marg_y.spines["top"].set_visible(False)
    ax_marg_y.spines["bottom"].set_visible(False)
    ax_marg_y.spines["left"].set_visible(False)
    ax_marg_y.spines["right"].set_visible(False)
    ax_marg_y.xaxis.set_ticks_position('none') 
    ax_marg_y.yaxis.set_ticks_position('none') 
    
    
    ax_marg_y.set_ylim((-0.05,1.05))
    ax_joint.set_ylim((-0.05,1.05))
    
    ax_marg_x.set_xlim((-0.15,3.05))
    ax_joint.set_xlim((-0.15,3.05))

    # Turn off tick labels on marginals
    plt.setp(ax_marg_x.get_xticklabels(), visible=False)
    plt.setp(ax_marg_x.get_yticklabels(), visible=False)
    plt.setp(ax_marg_y.get_xticklabels(), visible=False)
    plt.setp(ax_marg_y.get_yticklabels(), visible=False)
    
    
    xnewLabels = np.array([str(x) for x in xtags])
    ax_joint.set_xticks([np.log10(x) for x in xtags])
    ax_joint.set_xticklabels(xnewLabels)
    
    
    # Set labels on joint
    
    ax_joint.spines["right"].set_visible(False)
    ax_joint.spines["top"].set_visible(False)
#     ax_joint.spines['left'].set_linewidth(2)
#     ax_joint.spines['bottom'].set_linewidth(2)

#     ax_joint.xaxis.set_tick_params(width=2)
#     ax_joint.yaxis.set_tick_params(width=2)
    
    
    ######
    
    ax_joint.set_yticks([0, 0.5, 1])
    ax_joint.set_yticklabels(['0.0', '0.5', '1.0'])
    
    ax_joint.tick_params(labelsize=28, length=5)

    if not just_show:
        plt.savefig(plot_dir + '/' + event + '_' + dset + '.svg', bbox_inches='tight', transparent=True)
        plt.savefig(plot_dir + '/' + event + '_' + dset + '.pdf', bbox_inches='tight', transparent=True)
        plt.savefig(plot_dir + '/' + event + '_' + dset + '.png', bbox_inches='tight', transparent=True)
        
    else:
        plt.show()
        
    
def plot_histograms(hist_list, dset_name_list, ylab='ASE ranked by expression', plot_title = "Event $\Psi$ distributions",
                    plot_name = 'PSI_distributions', plot_dir = 'plots/', cmap='binary', fig_len=20, fig_height = 20, just_show=False,
                    ypos1=0, ypos2=0.6, ls=1, sk=1, tfs=16, name_size=20, nbins=25, hspace = None, wspace = None):
        
    if len(hist_list) != len(dset_name_list):
        raise Exception('Number of datasets and names do not match')
        return
    
    n_hist = len(hist_list)
    
    fig_width = int((fig_len/5) * n_hist)
    
    figsize(fig_width, fig_height)

    fig = gcf()
    
#     gs = GridSpec(1,4)
# gs.update(wspace=0.25, hspace=1.05)
    
#     if (hspace != '') and (wspace != ''):
#         gs.update(wspace=wspace, hspace=hspace)
    
    fig.suptitle(plot_title, fontsize=28)

    fig.text(ypos1, ypos2, ylab, va='center', rotation='vertical', fontsize=28)
    
    for i in range(n_hist-1):

        subplot(1, n_hist, i+1)
        plt.subplots_adjust(wspace=wspace, hspace=wspace)
        
        #.tick_params(labelsize=ls)

        ax = sns.heatmap(hist_list[i], cmap=cmap, xticklabels=nbins//2, yticklabels=False, 
                    cbar_kws = dict(use_gridspec=False,location="bottom"), cbar=False, fmt=".1f",
                    annot_kws={'size':28})
        plt.setp(ax.get_yticklabels(), visible=False)
        
        ax.set_xticks([0.5, nbins/2, nbins-0.5])
        ax.set_xticklabels(['0', '0.5', '1'])
        
        ax.tick_params(labelsize=28, length=5)
        
#         plt.setp(ax.get_yticklabels(), visible=False)

        plt.xlabel('$\hat{\Psi}$', fontsize=28)
        plt.title(dset_name_list[i], fontsize=name_size)
    
    
    subplot(1,n_hist,n_hist)
    
    #cbar_ax = fig.add_axes([.905, .3, .05, .3])

    ax = sns.heatmap(hist_list[-1], cmap=cmap, xticklabels=nbins//2, yticklabels=False, 
               cbar_kws = dict(use_gridspec=False,location="bottom",label='Relative frequency', shrink= sk))
#     ax.figure.axes[-1].yaxis.label.set_size(20)

    plt.xlabel('$\hat{\Psi}$', fontsize=28)
    plt.title(dset_name_list[-1], fontsize=name_size)
    
    ax.figure.axes[-1].yaxis.label.set_size(ls)
    ax.figure.axes[-1].xaxis.label.set_size(ls)
    
    plt.setp(ax.get_yticklabels(), visible=False)
        
    ax.set_xticks([0.5, nbins/2, nbins-0.5])
    ax.set_xticklabels(['0', '0.5', '1'])
        
    ax.tick_params(labelsize=28, length=5)
    
    ax.figure.axes[-1].tick_params(labelsize=28)
    
    #ax.figure.axes[-1].tick.label.set_size(20)
    #ax.figure.axes[-1].figsize=[10, 30]

    fig.tight_layout(rect=[0.05, 0.235, 1, 0.95])

    if not just_show:
#         plt.savefig(plot_dir + '/' + plot_name + '.svg')
        plt.savefig(plot_dir + '/' + plot_name + '.pdf', bbox_inches='tight', transparent=True)
        plt.savefig(plot_dir + '/' + plot_name + '.png', bbox_inches='tight', transparent=True)
        
    else:
        plt.show()

        
def plot_boxplot(cuentas_df, ylabel, dset_names, plot_name, xlocs = np.array([1, 2, 3, 4, 5]), 
                 ylocs=False, ytags=False, xlabel=False, save_dir='plots/analysis_3/', plot_line = False):
    
    fig = plt.figure(figsize=(6, 5))

    ax  = plt.subplot(1,1,1)
    
    if plot_line:
        
        plt.plot([0.75, np.max(xlocs)+0.25], [plot_line, plot_line], linestyle='--', c='darkorange', linewidth=1)
    
    bp = ax.boxplot([x for x in cuentas_df], showfliers=False)

    plt.setp(bp['boxes'], color='darkblue', linewidth=1)
    plt.setp(bp['whiskers'], color='darkblue', linewidth=1, linestyle='--')
    plt.setp(bp['caps'], color='black', linewidth=1)
    plt.setp(bp['medians'], color='darkred', linewidth=1)

    plt.ylabel(ylabel)
    if xlabel:
        plt.xlabel(xlabel)
    ax.figure.axes[-1].tick_params(labelsize=28, length=5)
    ax.figure.axes[-1].yaxis.label.set_size(28)
    ax.figure.axes[-1].xaxis.label.set_size(28)

    ax.xaxis.set_tick_params(width=1)
    ax.yaxis.set_tick_params(width=1)

    plt.xticks(rotation='vertical')

    #xlocs = np.array([1, 2, 3, 4, 5])
    #xlocs = np.array(range(len(dset_names)))+1
    xnewLabels = np.array(dset_names)
    plt.xticks(xlocs, xnewLabels)

    if ylocs and ytags:
        ylocs = np.array(ylocs)
        ynewLabels = np.array(ytags)
        plt.yticks(ylocs, ynewLabels)

    ax.tick_params(labelsize=28, length=5)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
#     ax.spines['left'].set_linewidth(3)
#     ax.spines['bottom'].set_linewidth(3)

#     ax.xaxis.set_tick_params(width=3)
#     ax.yaxis.set_tick_params(width=3)

    plt.savefig(save_dir + plot_name + '.svg', bbox_inches='tight', transparent=True)
    plt.savefig(save_dir + plot_name + '.pdf', bbox_inches='tight', transparent=True)
    plt.savefig(save_dir + plot_name + '.png', dpi=300, bbox_inches='tight', transparent=True)
    
    
    
def plot_splicing_space(PSI_tab, read_counts, mrna_counts, exon_list,
                        filter_exons=True, count_lim = 50):

    binary_exons = []
    mean_mrna = []
    mean_psi = []
    mean_psi_dist = []
    for e in exon_list:
        mean_mrna_counts = mrna_counts.loc[e.split('_')[0]].mean()
        if mean_mrna_counts <= 50:
            
            mean_psi_exon = PSI_tab.loc[e].mean()
            norm_factor = 1 + mean_psi_exon
            
            reads_per_molecule = reads_per_mrna(read_counts.loc[e], mrna_counts.loc[e.split('_')[0]], norm_factor) 
            if ((reads_per_molecule.quantile() <= 1000000) and (reads_per_molecule.quantile() >= 0.25)) or (not filter_exons):

                binary_exons.append((np.sum(PSI_tab.loc[e] == 1) + np.sum(PSI_tab.loc[e] == 0))/np.sum(~PSI_tab.loc[e].isnull()))
                mean_mrna.append(mrna_counts.loc[e.split('_')[0]].mean())
                mean_psi.append(mean_psi_exon)
                mean_psi_dist.append(0.5 - np.abs(0.5 - PSI_tab.loc[e].mean()))

    plt.figure()
    plt.title('mRNA counts and $\Psi$')
    sc = plt.scatter(mean_mrna, mean_psi_dist, c=np.array(binary_exons), cmap='bwr')
    plt.colorbar(sc, label='Proportion of binary observations')
    plt.ylabel('$mean \Psi$ (adjusted)')
    plt.xlabel('mRNA counts')

    
def plot_information_map(mrna_list, reads_list, binary_list, rpm_list, dset_name,
                         xlocs, xtags, ylocs, ytags, max_x = False, max_y = False,
                         stat = 'mean', color_by = 'binary', quant = 0.1, mrna_min = 15, set_sr_lim = False,
                         save_name = '', plot_dir = 'plots/analysis_3/'):
    
    fig = plt.figure(figsize=(7, 6))
    ax  = plt.subplot(1,1,1)
    
    if color_by == 'binary':
        cmap = 'bwr'
    elif color_by == 'psi':
        colors1 = plt.cm.viridis_r(np.linspace(0., 1, 128))
        colors2 = plt.cm.viridis(np.linspace(0, 1, 128))

        # combine them and build a new colormap
        colors = np.vstack((colors1, colors2))
        cmap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
    else:
        raise Exception('invalid color map')

    sc = plt.scatter(np.log10(np.array(mrna_list)+1), np.log10(np.array(reads_list)+1), 
                     c=binary_list, alpha=0.75, cmap=cmap, s=75, edgecolors='none', vmin=0, vmax=1)
    
    if (not max_x) and (not max_y):
        max_x = np.int(np.ceil(np.log10(np.array(mrna_list)+1).max()))
        max_y = np.int(np.ceil(np.log10(np.array(reads_list)+1).max()))

        max_x = np.max((max_y, max_x))
        max_y = max_x
    
    if set_sr_lim:
        min_sr = set_sr_lim
    else:
        if stat == 'mean':
            min_sr = rpm_list.mean(axis=1).quantile(0.1)
        elif stat == 'median':
            min_sr = rpm_list.median(axis=1).quantile(0.1)
        else:
            raise Exception('invalid stat parameter')
    
    #ax.plot(np.array([0,np.min((max_x, max_y))]), np.array([0,np.min((max_x, max_y))]), c='darkred')
    #ax.plot(np.array([0,np.min((max_x, max_y))]), np.array([0,np.min((max_x, max_y))]), c='darkred')
    
#     ax.fill_between(np.log10(np.arange(0,10**max_x, 0.1) + 1),
#                          np.log10(np.arange(0,10**max_x, 0.1)*1.95 + 1), 
#                          np.log10(np.arange(0,10**max_x, 0.1)*1.05 + 1), color='navy', zorder=3, alpha=0.2, linewidth=0.0)
    
    ax.plot(np.log10(np.arange(0,10**max_x, 0.1) + 1), np.log10(np.arange(0,10**max_x, 0.1)*rpm_list.median(axis=1).median() + 1), '--', c='darkred')
    
    
    ax.fill_between(np.log10(np.arange(0,10**max_x, 0.1) + 1),
                         np.log10(np.arange(0,10**max_x, 0.1)*min_sr + 1), 
                         np.zeros(10**(max_x+1)), color='darkred', zorder=3, alpha=0.3, linewidth=0.0)
    
    ax.fill_between([0, np.log10(mrna_min+1)],
                     [np.log10(1.5*(10**max_y)),np.log10(1.5*(10**max_y))], [0,0], color='darkred', zorder=3, alpha=0.3, linewidth=0.0)
    #ax.plot([0,3], [0,3], c='darkred')
    plt.title(dset_name, fontsize =28)
    plt.ylabel('average SJ reads', fontsize=28)
    plt.xlabel('average mRNA counts', fontsize=28)
    cb = plt.colorbar(sc)
    cb.set_label(label='Proportion of binary $\Psi$', size=28)
    cb.ax.tick_params(labelsize=28, length=5)

    cb.outline.set_visible(False)

    ax.tick_params(labelsize=28, length=5)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
#     ax.spines['left'].set_linewidth(3)
#     ax.spines['bottom'].set_linewidth(3)

#     ax.xaxis.set_tick_params(width=2)
#     ax.yaxis.set_tick_params(width=2)
    
    #xlocs = [0] + [int(np.log10(10**i+1)+0.1) for i in range(max_x+1)]
    #xtags = [str(10**x) for x in xlocs]

    xlocs = np.array(xlocs)
    xnewLabels = np.array(xtags)
    plt.xticks(xlocs, xnewLabels)
    
    #ylocs = [0] + [int(np.log10(10**i+1)+0.1) for i in range(max_y+1)]
    #ytags = [str(10**x) for x in ylocs]

    ylocs = np.array(ylocs)
    ynewLabels = np.array(ytags)
    plt.yticks(ylocs, ynewLabels)
    
    plt.savefig(plot_dir + save_name + '_information_map.svg', bbox_inches='tight', transparent=True)
    plt.savefig(plot_dir + save_name + '_information_map.pdf', bbox_inches='tight', transparent=True)
    plt.savefig(plot_dir + save_name + '_information_map.png', dpi=300, bbox_inches='tight', transparent=True)

    
    
def plot_information_map2(mrna_list, reads_list, binary_list, dset_name, # removed rpm_list; unnecessary
                         xlocs, xtags, ylocs, ytags, max_x = False, max_y = False,
                         stat = 'mean', color_by = 'binary', quant = 0.1, mrna_min = 10, read_min=10, #set_sr_lim = False,
                         save_name = '', plot_dir = 'plots/analysis_5/'):
    
    fig = plt.figure(figsize=(7, 5))
    ax  = plt.subplot(1,1,1)
    
    if color_by == 'binary':
        cmap = 'bwr'
    elif color_by == 'psi':
        colors1 = plt.cm.viridis_r(np.linspace(0., 1, 128))
        colors2 = plt.cm.viridis(np.linspace(0, 1, 128))

        # combine them and build a new colormap
        colors = np.vstack((colors1, colors2))
        cmap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
    else:
        raise Exception('invalid color map')

    sc = plt.scatter(np.log10(np.array(mrna_list)+1), np.log10(np.array(reads_list)+1), 
                     c=binary_list, alpha=0.75, cmap=cmap, s=75, edgecolors='none', vmin=0, vmax=1)
    
    if (not max_x) and (not max_y):
        max_x = np.int(np.ceil(np.log10(np.array(mrna_list)+1).max()))
        max_y = np.int(np.ceil(np.log10(np.array(reads_list)+1).max()))

        max_x = np.max((max_y, max_x))
        max_y = max_x
    
#     if set_sr_lim:
#         min_sr = set_sr_lim
#     else:
#         if stat == 'mean':
#             min_sr = rpm_list.mean(axis=1).quantile(0.1)
#         elif stat == 'median':
#             min_sr = rpm_list.median(axis=1).quantile(0.1)
#         else:
#             raise Exception('invalid stat parameter')
    
    #ax.plot(np.array([0,np.min((max_x, max_y))]), np.array([0,np.min((max_x, max_y))]), c='darkred')
    #ax.plot(np.array([0,np.min((max_x, max_y))]), np.array([0,np.min((max_x, max_y))]), c='darkred')
    
#     ax.fill_between(np.log10(np.arange(0,10**max_x, 0.1) + 1),
#                          np.log10(np.arange(0,10**max_x, 0.1)*1.95 + 1), 
#                          np.log10(np.arange(0,10**max_x, 0.1)*1.05 + 1), color='navy', zorder=3, alpha=0.2, linewidth=0.0)
    
    #ax.plot(np.log10(np.arange(0,10**max_x, 0.1) + 1), np.log10(np.arange(0,10**max_x, 0.1)*rpm_list.median(axis=1).median() + 1), '--', c='darkred')
    
#     ax.plot([np.log10(mrna_min+1), np.log10(mrna_min+1)], [0, np.log10(1.1*(10**max_y))], '--', c='black')
#     ax.plot([0, np.log10(10**max_x + 1)], [np.log10(read_min+1), np.log10(read_min+1)], '--', c='black')
    
    
#     ax.fill_between([0,np.log10(10**max_x + 1)], [np.log10(read_min+1), np.log10(read_min+1)], [0,0], 
#                     color='gray', zorder=3, alpha=0.3, linewidth=0.0)
    
#     ax.fill_between([0, np.log10(mrna_min+1)],
#                      [np.log10(1.1*(10**max_y)),np.log10(1.1*(10**max_y))], [0,0], color='gray', zorder=3, alpha=0.3, linewidth=0.0)
    
#     #ax.plot([0,3], [0,3], c='darkred')
    plt.title(dset_name, fontsize =28)
    plt.ylabel('average SJ reads', fontsize=28)
    plt.xlabel('average mRNA counts', fontsize=28)
    cb = plt.colorbar(sc)
    cb.set_label(label='Proportion of binary $\Psi$', size=28)
    cb.ax.tick_params(labelsize=28, length=5)

    cb.outline.set_visible(False)

    ax.tick_params(labelsize=28, length=5)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
#     ax.spines['left'].set_linewidth(3)
#     ax.spines['bottom'].set_linewidth(3)

#     ax.xaxis.set_tick_params(width=2)
#     ax.yaxis.set_tick_params(width=2)
    
    #xlocs = [0] + [int(np.log10(10**i+1)+0.1) for i in range(max_x+1)]
    #xtags = [str(10**x) for x in xlocs]

    xlocs = np.array(xlocs)
    xnewLabels = np.array(xtags)
    plt.xticks(xlocs, xnewLabels)
    
    #ylocs = [0] + [int(np.log10(10**i+1)+0.1) for i in range(max_y+1)]
    #ytags = [str(10**x) for x in ylocs]

    ylocs = np.array(ylocs)
    ynewLabels = np.array(ytags)
    plt.yticks(ylocs, ynewLabels)
    
    plt.savefig(plot_dir + save_name + '_information_map.svg', bbox_inches='tight', transparent=True)
    plt.savefig(plot_dir + save_name + '_information_map.pdf', bbox_inches='tight', transparent=True)
    plt.savefig(plot_dir + save_name + '_information_map.png', dpi=300, bbox_inches='tight', transparent=True)



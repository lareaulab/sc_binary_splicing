from tqdm import tqdm_notebook as tqdm
current_palette = sns.color_palette('dark')
from sklearn.metrics import auc
from scipy.stats.mstats import mannwhitneyu
from scipy.stats import wilcoxon


def get_bimodality(PSI_tab, x):
    extreme_cases = ((PSI_tab <= x) | (PSI_tab >= (1-x))).sum(axis=1)# + (PSI_tab >= (1-x)).sum(axis=1))
    total_cases = (~PSI_tab.isna()).sum(axis=1)
    bimodality = extreme_cases / total_cases
    return bimodality

def get_bimodality_curves(PSI_tab, steps):
    bimodality_curves = pd.DataFrame()
    for x in np.arange(0, 0.5+steps, steps):
#         x = np.arange(0, 0.51, 0.01)[i]
        
        bimodality_x = get_bimodality(PSI_tab, x)
        bimodality_curves[str(round(x, 3))] = bimodality_x
        
    bimodality_curves.index = PSI_tab.index
    return bimodality_curves

def get_curves_dataset(PSI_tab, subpop, mrna_counts, mrna_per_event, read_counts, 
                       coverage_tab, psi_min = 0.2, mrna_min = 10, reads_min = 0, cell_min = 0.5, filter_cj = True,
                      get_all = True, steps = 0.01):
    
    int_ = PSI_tab[subpop].loc[np.abs(0.5 - PSI_tab[subpop].mean(axis=1)) <= (0.5-psi_min)]
    PSI_int = int_.loc[int_.isna().mean(axis=1) <= (1-cell_min)]

    if get_all:
        bimodality_curves_all = get_bimodality_curves(PSI_int, steps)
        
    #########
    
    PSI_filtered = process_subpop(subpop, PSI_tab, mrna_counts, 
                                        mrna_per_event, read_counts, coverage_tab['SJ_coverage'], 
                                       psi_min, mrna_min, reads_min = reads_min, cell_min=cell_min, filter_cj = filter_cj)
    
#     assert np.all([x in PSI_int.index for x in PSI_filtered[0].index])
    
    bimodality_curves_selected =get_bimodality_curves(PSI_tab.loc[PSI_filtered[0].index], steps)
    bimodality_curves_filtered =get_bimodality_curves(PSI_filtered[0], steps)

    if get_all:
        return bimodality_curves_all, bimodality_curves_selected, bimodality_curves_filtered
    
    else:
        return bimodality_curves_selected, bimodality_curves_filtered
        
        
        
def plot_bimodality_curves(curva_set, color_list, nameset, 
                            title, alpha_list = '', plot_individual_events = False):
    
    if len(alpha_list) == 0:
        alpha_list = [0.5]*len(curva_set)
    
    for i in range(len(curva_set)):
        
        bimodality_curve = curva_set[i]
        col = color_list[i]
        alph = alpha_list[i]
        x_lims = [float(x) for x in bimodality_curve.columns]
        
        
        for exon in bimodality_curve.index:
            curva = list(bimodality_curve.loc[exon])
            plt.plot([0]+x_lims, [0]+list(curva), c=col, alpha=alph,)
            
    plt.title(title)
    plt.xlabel('$\Psi$ cutoff')
    plt.ylabel('proportion of cells below cutoff')

    colors = color_list
    lines = [Line2D([0], [0], color=c, linewidth=3) for c in colors]
    labels = nameset
    plt.legend(lines, labels, loc='lower right', frameon=False)

    plt.show()
    
# def plot_bimodality_box(curva_set, color_list, nameset, 
#                             title, alpha_list = '', plot_individual_events = False):
    
#     if len(alpha_list) == 0:
#         alpha_list = [0.5]*len(curva_set)
    
#     for i in range(len(curva_set)):
        
#         bimodality_curve = curva_set[i]
#         col = color_list[i]
#         alph = alpha_list[i]
#         x_lims = [float(x) for x in bimodality_curve.columns]
        
        
#         for exon in bimodality_curve.index:
#             curva = list(bimodality_curve.loc[exon])
#             plt.plot([0]+x_lims, [0]+list(curva), c=col, alpha=alph,)
            
#     plt.title(title)
#     plt.xlabel('$\Psi$ cutoff')
#     plt.ylabel('proportion of cells below cutoff')

#     colors = color_list
#     lines = [Line2D([0], [0], color=c, linewidth=3) for c in colors]
#     labels = nameset
#     plt.legend(lines, labels, loc='lower right', frameon=False)

#     plt.show()
    
def process_bimodality(PSI_tab, subpop, mrna_counts, mrna_per_event, read_counts, coverage_tab, dset, psi_min = 0.2, 
                       mrna_min = 10, reads_min = 10, cell_min = 0.5, steps = 0.05):
    
    curves_all, mrna_selected, mrna_filtered = get_curves_dataset(PSI_tab, subpop, mrna_counts, mrna_per_event, read_counts, 
                       coverage_tab, psi_min = psi_min, mrna_min = mrna_min, reads_min = 0, cell_min = cell_min, steps=steps)

    
    read_selected, read_filtered = get_curves_dataset(PSI_tab, subpop, mrna_counts, mrna_per_event, read_counts, 
                       coverage_tab, psi_min = psi_min, mrna_min = 0, reads_min = reads_min, 
                                                      cell_min = cell_min, steps=steps, get_all=False)
    
    
    mrna_only_selected, mrna_only_filtered = get_curves_dataset(PSI_tab, subpop, mrna_counts, mrna_per_event, read_counts, 
                       coverage_tab, psi_min = psi_min, mrna_min = mrna_min, reads_min = 0, cell_min = cell_min, steps=steps,
                                                                     filter_cj = False, get_all = False)
    
    
    plot_bimodality_boxplots(curves_all, read_filtered, mrna_only_filtered, mrna_filtered, dset)
    
    plot_bimodality_difference(curves_all, read_filtered, mrna_only_filtered, mrna_filtered, dset)

    
    run_AUC_tests(curves_all, read_filtered, mrna_only_filtered, mrna_filtered)
    
    
    
def get_auc(bimodality_curves):
    
    bimodal_limits = [float(x) for x in bimodality_curves.columns]
    
    auc_list = []
    for exon in bimodality_curves.index:
        auc_list.append(auc(bimodal_limits, bimodality_curves.loc[exon]))
        
    return auc_list


def run_AUC_tests(curves_all, read_filtered, mrna_only_filtered, mrna_filtered):
    
    auc_all = get_auc(curves_all)
    
    auc_read_selected = get_auc(curves_all.loc[read_filtered.index])
    auc_read_filtered = get_auc(read_filtered)
    
    auc_mrna_only_selected = get_auc(curves_all.loc[mrna_only_filtered.index])
    auc_mrna_only_filtered = get_auc(mrna_only_filtered)
    
    auc_mrna_selected = get_auc(curves_all.loc[mrna_filtered.index])
    auc_mrna_filtered = get_auc(mrna_filtered)
    
    pval_read_selected = mannwhitneyu(auc_all, auc_read_filtered)[1]
    pval_mrna_only_selected = mannwhitneyu(auc_all, auc_mrna_only_filtered)[1]
    pval_mrna_selected = mannwhitneyu(auc_all, auc_mrna_filtered)[1]
    
    pval_read_filtered = mannwhitneyu(auc_read_selected, auc_read_filtered)[1]
    pval_mrna_only_filtered = mannwhitneyu(auc_mrna_only_selected, auc_mrna_only_filtered)[1]
    pval_mrna_filtered = mannwhitneyu(auc_mrna_selected, auc_mrna_filtered)[1]
    
    pvals_adj = multipletests([pval_read_selected, pval_mrna_only_selected, pval_mrna_selected,
                               pval_read_filtered, pval_mrna_only_filtered, pval_mrna_filtered], method='fdr_bh')
    
    
    print('Median AUC in all exons: ' + str(median(auc_all)))
    print('')
    print('Median AUC in read selected exons: ' + str(median(auc_read_filtered)))
    if pvals_adj[0][0]:
        print('Change is significant with respect to all exons (adj pval = ' + str(pvals_adj[1][0]) + ')')
    else:
        print('Change is not significant with respect to all exons (adj pval = ' + str(pvals_adj[1][0]) + ')')
        
    if pvals_adj[0][3]:
        print('Filtering observations reduces AUC significantly (adj pval = ' + str(pvals_adj[1][3]) + ')')
    else:
        print('Change is not significant with respect to all exons (adj pval = ' + str(pvals_adj[1][3]) + ')')
        
    print('')
        
        
        
        
    print('Median AUC in mRNA only selected exons: ' + str(median(auc_mrna_only_filtered)))
    if pvals_adj[0][1]:
        print('Change is significant with respect to all exons (adj pval = ' + str(pvals_adj[1][1]) + ')')
    else:
        print('Filtering observations does not reduce AUC significantly (adj pval = ' + str(pvals_adj[1][1]) + ')')
        
    if pvals_adj[0][4]:
        print('Filtering observations reduces AUC significantly (adj pval = ' + str(pvals_adj[1][4]) + ')')
    else:
        print('Filtering observations does not reduce AUC significantly (adj pval = ' + str(pvals_adj[1][4]) + ')')
        
    print('')
    
    
    print('Median AUC in mRNA selected exons: ' + str(median(auc_mrna_filtered)))
    if pvals_adj[0][2]:
        print('Change is significant with respect to all exons (adj pval = ' + str(pvals_adj[1][2]) + ')')
    else:
        print('Change is not significant with respect to all exons (adj pval = ' + str(pvals_adj[1][2]) + ')')
        
    if pvals_adj[0][5]:
        print('Filtering observations reduces AUC significantly (adj pval = ' + str(pvals_adj[1][5]) + ')')
    else:
        print('Filtering observations does not reduce AUC significantly (adj pval = ' + str(pvals_adj[1][5]) + ')')
        
    print('')
        
    
#     return pvals_adj[1]

####################################################################

def run_AUC_tests(curves_all, read_filtered, mrna_only_filtered, mrna_filtered):
    
    auc_all = get_auc(curves_all)
    
    auc_read_selected = get_auc(curves_all.loc[read_filtered.index])
    auc_read_filtered = get_auc(read_filtered)
    
    auc_mrna_only_selected = get_auc(curves_all.loc[mrna_only_filtered.index])
    auc_mrna_only_filtered = get_auc(mrna_only_filtered)
    
    auc_mrna_selected = get_auc(curves_all.loc[mrna_filtered.index])
    auc_mrna_filtered = get_auc(mrna_filtered)
    
    pval_read_selected = mannwhitneyu(auc_all, auc_read_filtered)[1]
    pval_mrna_only_selected = mannwhitneyu(auc_all, auc_mrna_only_filtered)[1]
    pval_mrna_selected = mannwhitneyu(auc_all, auc_mrna_filtered)[1]
    
    pval_read_filtered = mannwhitneyu(auc_read_selected, auc_read_filtered)[1]
    pval_mrna_only_filtered = mannwhitneyu(auc_mrna_only_selected, auc_mrna_only_filtered)[1]
    pval_mrna_filtered = mannwhitneyu(auc_mrna_selected, auc_mrna_filtered)[1]
    
    pvals_adj = multipletests([pval_read_selected, pval_mrna_only_selected, pval_mrna_selected,
                               pval_read_filtered, pval_mrna_only_filtered, pval_mrna_filtered], method='fdr_bh')
    
    
    print('Median AUC in all exons: ' + str(median(auc_all)))
    print('')
    print('Median AUC in read selected exons: ' + str(median(auc_read_filtered)))
    if pvals_adj[0][0]:
        print('Change is significant with respect to all exons (adj pval = ' + str(pvals_adj[1][0]) + ')')
    else:
        print('Change is not significant with respect to all exons (adj pval = ' + str(pvals_adj[1][0]) + ')')
        
    if pvals_adj[0][3]:
        print('Filtering observations reduces AUC significantly (adj pval = ' + str(pvals_adj[1][3]) + ')')
    else:
        print('Change is not significant with respect to all exons (adj pval = ' + str(pvals_adj[1][3]) + ')')
        
    print('')
        
        
        
        
    print('Median AUC in mRNA only selected exons: ' + str(median(auc_mrna_only_filtered)))
    if pvals_adj[0][1]:
        print('Change is significant with respect to all exons (adj pval = ' + str(pvals_adj[1][1]) + ')')
    else:
        print('Filtering observations does not reduce AUC significantly (adj pval = ' + str(pvals_adj[1][1]) + ')')
        
    if pvals_adj[0][4]:
        print('Filtering observations reduces AUC significantly (adj pval = ' + str(pvals_adj[1][4]) + ')')
    else:
        print('Filtering observations does not reduce AUC significantly (adj pval = ' + str(pvals_adj[1][4]) + ')')
        
    print('')
    
    
    print('Median AUC in mRNA selected exons: ' + str(median(auc_mrna_filtered)))
    if pvals_adj[0][2]:
        print('Change is significant with respect to all exons (adj pval = ' + str(pvals_adj[1][2]) + ')')
    else:
        print('Change is not significant with respect to all exons (adj pval = ' + str(pvals_adj[1][2]) + ')')
        
    if pvals_adj[0][5]:
        print('Filtering observations reduces AUC significantly (adj pval = ' + str(pvals_adj[1][5]) + ')')
    else:
        print('Filtering observations does not reduce AUC significantly (adj pval = ' + str(pvals_adj[1][5]) + ')')
        
    print('')
        


def process_bimodality_dset(PSI_tab, subpop_list, subpop_names, mrna_counts, mrna_per_event, read_counts, coverage_tab, 
                            dset, psi_min = 0.2, mrna_min = 10, mrna_read_min=0, reads_min = 10, cell_min = 0.5, steps = 0.25):
    
    subpop_all = []
    subpop_mrna_filtered = []
    subpop_read_filtered = []
    
    for i in range(len(subpop_list)):
        subpop = subpop_list[i]
        curves_all, mrna_selected, mrna_filtered = get_curves_dataset(PSI_tab, subpop, mrna_counts, mrna_per_event, read_counts, 
                           coverage_tab, psi_min = psi_min, mrna_min = mrna_min, reads_min = mrna_read_min, 
                                                                      cell_min = cell_min, steps=steps)


        read_selected, read_filtered = get_curves_dataset(PSI_tab, subpop, mrna_counts, mrna_per_event, read_counts, 
                           coverage_tab, psi_min = psi_min, mrna_min = 0, reads_min = reads_min, 
                                                          cell_min = cell_min, steps=steps, get_all=False)
        
#         mrna_only_selected, mrna_only_filtered = get_curves_dataset(PSI_tab, subpop, mrna_counts, mrna_per_event, read_counts, 
#                            coverage_tab, psi_min = psi_min, mrna_min = mrna_min, reads_min = 0, cell_min = cell_min, steps=steps,
#                                                                          filter_cj = False, get_all = False)

        print(curves_all.columns[1])
        subpop_all.append(list(curves_all[curves_all.columns[1]]))
        subpop_mrna_filtered.append(list(mrna_filtered[mrna_filtered.columns[1]]))
        subpop_read_filtered.append(list(read_filtered[read_filtered.columns[1]]))


    
    plot_bimodality_boxplots_dset(subpop_all, subpop_read_filtered, subpop_mrna_filtered, subpop_list, subpop_names, dset)
    
    
#     run_AUC_tests(curves_all, read_filtered, mrna_only_filtered, mrna_filtered)
    
    
def plot_bimodality_boxplots_dset(curves_all, read_selected, mrna_selected, subpop_list, subpop_names, dset):

#     fig = plt.figure(figsize=(((9-((4.5/4)*(4-len(curves_all))))/4)*len(curves_all),6))

    fig = plt.figure(figsize=( ( ((4.2*9)/(5*4))*len(curves_all) + (((len(curves_all)-1))) ,4.5)))
    ax  = plt.subplot(1,1,1)

    data_a = curves_all
    data_b = read_selected
    data_c = mrna_selected

    # ticks = ['A', 'B', 'C']

    # def set_box_color(bp, color):
    # #     plt.setp(bp['face_color'], color=color)
    # #     plt.setp(bp['whiskers'], color=color)
    # #     plt.setp(bp['caps'], color=color)
    #     plt.setp(bp['medians'], color=darkred)

    # plt.figure()

    bpl = ax.boxplot(data_a, positions=np.array(range(len(data_a)))*5.0-1.15, sym='', widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor=current_palette[4], color=current_palette[4]))
    bpr = ax.boxplot(data_b, positions=np.array(range(len(data_b)))*5.0, sym='', widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='darkred', color='darkred'))
    bpm = ax.boxplot(data_c, positions=np.array(range(len(data_c)))*5.0+1.15, sym='', widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='navy', color='navy'))


    [x.set_facecolor(current_palette[4]) for x in bpl['boxes']]
    [x.set_facecolor(current_palette[2]) for x in bpr['boxes']]
#     [x.set_facecolor(current_palette[1]) for x in bpmo['boxes']]
    [x.set_facecolor(current_palette[0]) for x in bpm['boxes']]
    # plt.setp(bpl['medians'], color='darkorange')
    # plt.setp(bpr['medians'], color='darkorange')

    # for patch, color in zip(bplot1['boxes'], colors):
    #     patch.set_facecolor(color)

    # set_box_color(bpl, 'lightblue') # colors are from http://colorbrewer2.org/
    # set_box_color(bpr, 'lightred')

    # draw temporary red and blue lines and use them to create a legend
    
    ax.tick_params(labelsize=28, length=5)###
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    
    plt.xlabel('Bimodality limit', fontsize=28)
    plt.ylabel('Percent bimodal cells', fontsize=28)
    plt.title(dset, fontsize=28)

    plt.xticks(np.array(range(len(data_b)))*5.0, subpop_names)
#     plt.xticks(range(0, len(curves_all.columns[:-1]) * 7, 14), curves_all.columns[:-1][::2])
    # plt.xlim(-2, len(ticks)*2)
    plt.tight_layout()
#     plt.xlim(-3, len(data_b)*7.0)
    plt.ylim(-0.05, 1.05)
    
    plt.yticks([0, 0.5, 1], ['0.0', '0.5', '1.0'])
    
    
    ###
#     plt.tick_params(labelsize=34, length=5)
    ###
    
    plot_name = '_'.join(dset.split())
    
    
    plt.savefig('plots_review/figure3/bimodality_plots_' + plot_name + '_subpop.svg', bbox_inches='tight', transparent=True)
    plt.savefig('plots_review/figure3/bimodality_plots_' + plot_name + '_subpop.pdf', bbox_inches='tight', transparent=True)
    plt.savefig('plots_review/figure3/bimodality_plots_' + plot_name + '_subpop.png', dpi=300, bbox_inches='tight', transparent=True)
    
    
    
    figsize(4,3)
    fig = plt.figure()
    ax  = plt.subplot(1,1,1)

    plt.plot([], c=current_palette[4], label='all exons', linewidth=5)
    plt.plot([], c=current_palette[2], label='reads filter', linewidth=5)
#     plt.plot([], c=current_palette[1], label='mRNA only filter', linewidth=5)
    plt.plot([], c=current_palette[0], label='mRNA filter', linewidth=5)
    plt.legend(frameon=False, fontsize=28)

    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)
    ax.xaxis.set_ticks_position('none') 
    ax.yaxis.set_ticks_position('none')
    plt.legend(fontsize=28, frameon=False)

    plt.savefig('plots_review/figure3/bimodality_plots_labels.svg', bbox_inches='tight', transparent=True)
    plt.savefig('plots_review/figure3/bimodality_plots_labels.pdf', bbox_inches='tight', transparent=True)
    plt.savefig('plots_review/figure3/bimodality_plots_labels.png', bbox_inches='tight', transparent=True)
    plt.show()
    
    
    
    
###############################################################    
# figsize(12, 6)

# figsize(12, 6)

def plot_bimodality_boxplots(curves_all, read_selected, mrna_only_selected, mrna_selected, dset):

    fig = plt.figure(figsize=(9,6))
    ax  = plt.subplot(1,1,1)

    data_a = [curves_all[x] for x in curves_all.columns[:-1]]
    data_b = [read_selected[x] for x in read_selected.columns[:-1]]
    data_c = [mrna_only_selected[x] for x in mrna_only_selected.columns[:-1]]
    data_d = [mrna_selected[x] for x in mrna_selected.columns[:-1]]

    # ticks = ['A', 'B', 'C']

    # def set_box_color(bp, color):
    # #     plt.setp(bp['face_color'], color=color)
    # #     plt.setp(bp['whiskers'], color=color)
    # #     plt.setp(bp['caps'], color=color)
    #     plt.setp(bp['medians'], color=darkred)

    # plt.figure()

    bpl = ax.boxplot(data_a, positions=np.array(range(len(data_a)))*7.0-1.5, sym='', widths=0.6, patch_artist=True)
    bpr = ax.boxplot(data_b, positions=np.array(range(len(data_b)))*7.0-0.5, sym='', widths=0.6, patch_artist=True)
    bpmo = ax.boxplot(data_c, positions=np.array(range(len(data_c)))*7.0+0.5, sym='', widths=0.6, patch_artist=True)
    bpm = ax.boxplot(data_d, positions=np.array(range(len(data_d)))*7.0+1.5, sym='', widths=0.6, patch_artist=True)


    [x.set_facecolor(current_palette[4]) for x in bpl['boxes']]
    [x.set_facecolor(current_palette[2]) for x in bpr['boxes']]
    [x.set_facecolor(current_palette[1]) for x in bpmo['boxes']]
    [x.set_facecolor(current_palette[0]) for x in bpm['boxes']]
    # plt.setp(bpl['medians'], color='darkorange')
    # plt.setp(bpr['medians'], color='darkorange')

    # for patch, color in zip(bplot1['boxes'], colors):
    #     patch.set_facecolor(color)

    # set_box_color(bpl, 'lightblue') # colors are from http://colorbrewer2.org/
    # set_box_color(bpr, 'lightred')

    # draw temporary red and blue lines and use them to create a legend
    
    ax.tick_params(labelsize=28, length=5)###
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    
    plt.xlabel('Bimodality limit', fontsize=28)
    plt.ylabel('Percent bimodal cells', fontsize=28)
    plt.title(dset, fontsize=28)

#     plt.xticks(range(0, len(curves_all.columns[:-1]) * 7, 7), curves_all.columns[:-1])
    plt.xticks(range(0, len(curves_all.columns[:-1]) * 7, 14), curves_all.columns[:-1][::2])
    # plt.xlim(-2, len(ticks)*2)
    plt.tight_layout()
    plt.xlim(-3, len(curves_all.columns[:-1])*7-3)
    plt.ylim(-0.05, 1.05)
    
    plt.yticks([0, 0.5, 1], ['0.0', '0.5', '1.0'])
    
    
    ###
#     plt.tick_params(labelsize=34, length=5)
    ###
    
    plot_name = '_'.join(dset.split())
    
    
    plt.savefig('plots_review/figure3/bimodality_plots_' + plot_name + '.svg', bbox_inches='tight')
    plt.savefig('plots_review/figure3/bimodality_plots_' + plot_name + '.pdf', bbox_inches='tight')
    plt.savefig('plots_review/figure3/bimodality_plots_' + plot_name + '.png', dpi=300, bbox_inches='tight')
    
    
    
    figsize(4,3)
    fig = plt.figure()
    ax  = plt.subplot(1,1,1)

    plt.plot([], c=current_palette[4], label='all exons', linewidth=5)
    plt.plot([], c=current_palette[2], label='reads filter', linewidth=5)
    plt.plot([], c=current_palette[1], label='mRNA only filter', linewidth=5)
    plt.plot([], c=current_palette[0], label='mRNA filter', linewidth=5)
    plt.legend(frameon=False, fontsize=28)

    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)
    ax.xaxis.set_ticks_position('none') 
    ax.yaxis.set_ticks_position('none')
    plt.legend(fontsize=28, frameon=False)

    plt.savefig('plots_review/figure3/bimodality_plots_labels.svg', bbox_inches='tight')
    plt.savefig('plots_review/figure3/bimodality_plots_labels.pdf', bbox_inches='tight')
    plt.savefig('plots_review/figure3/bimodality_plots_labels.png', bbox_inches='tight')
    plt.show()
    
    
def plot_bimodality_difference(curves_all, read_filtered, mrna_only_filtered, mrna_filtered, dset):

    fig = plt.figure(figsize=(9,6))
    ax  = plt.subplot(1,1,1)

    # data_a = [curves_all[x] for x in curves_all.columns[:-1]]
    data_a = [read_filtered[x] - curves_all.loc[read_filtered.index, x] for x in read_filtered.columns[:-1]]
    data_b = [mrna_only_filtered[x] - curves_all.loc[mrna_only_filtered.index, x] for x in mrna_only_filtered.columns[:-1]]
    data_c = [mrna_filtered[x] - curves_all.loc[mrna_filtered.index, x] for x in mrna_filtered.columns[:-1]]
    # ticks = ['A', 'B', 'C']


    bpr = ax.boxplot(data_a, positions=np.array(range(len(data_a)))*7.0-1.5, sym='', widths=0.9, patch_artist=True)
    bpmo = ax.boxplot(data_b, positions=np.array(range(len(data_b)))*7.0, sym='', widths=0.9, patch_artist=True)
    bpm = ax.boxplot(data_c, positions=np.array(range(len(data_c)))*7.0+1.5, sym='', widths=0.9, patch_artist=True)


    # [x.set_facecolor(current_palette[0]) for x in bpl['boxes']]
    [x.set_facecolor(current_palette[2]) for x in bpr['boxes']]
    [x.set_facecolor(current_palette[1]) for x in bpmo['boxes']]
    [x.set_facecolor(current_palette[0]) for x in bpm['boxes']]

    # plt.setp(bpl['medians'], color='darkorange')
    # plt.setp(bpr['medians'], color='darkorange')

    # for patch, color in zip(bplot1['boxes'], colors):
    #     patch.set_facecolor(color)

    # set_box_color(bpl, 'lightblue') # colors are from http://colorbrewer2.org/
    # set_box_color(bpr, 'lightred')

    # draw temporary red and blue lines and use them to create a legend
#     plt.plot([], c=current_palette[4], label='all exons')
    
#     plt.legend(frameon=False, fontsize=24)
    ax.tick_params(labelsize=28, length=5)##
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    
    plt.xlabel('Bimodality limit', fontsize=28)
    plt.ylabel('Bimodality change', fontsize=28)
    plt.title(dset, fontsize=28)

    plt.xticks(range(0, len(curves_all.columns[:-1]) * 7, 14), curves_all.columns[:-1][::2])
#     plt.xticks(range(0, len(curves_all.columns[:-1]) * 7, 7), curves_all.columns[:-1])
    
    plt.tight_layout()
    plt.xlim(-3, len(curves_all.columns[:-1])*7-3)
#     plt.ylim(-0.05, 1.05)
    
    plot_name = '_'.join(dset.split())
    
    
    plt.savefig('plots_review/figure3/bimodality_reduction_' + plot_name + '.svg', bbox_inches='tight')
    plt.savefig('plots_review/figure3/bimodality_reduction_' + plot_name + '.pdf', bbox_inches='tight')
    plt.savefig('plots_review/figure3/bimodality_reduction_' + plot_name + '.png', dpi=300, bbox_inches='tight')
    
    
    figsize(4,3)
    fig = plt.figure()
    ax  = plt.subplot(1,1,1)

    plt.plot([], c=current_palette[2], label='reads filter', linewidth=5)
    plt.plot([], c=current_palette[1], label='mRNA only filter', linewidth=5)
    plt.plot([], c=current_palette[0], label='mRNA filter', linewidth=5)
    plt.legend(frameon=False, fontsize=28)

    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)
    ax.xaxis.set_ticks_position('none') 
    ax.yaxis.set_ticks_position('none')
    plt.legend(fontsize=28, frameon=False)

    plt.savefig('plots_review/figure3/bimodality_reduction_labels.svg', bbox_inches='tight')
    plt.savefig('plots_review/figure3/bimodality_reduction_labels.pdf', bbox_inches='tight')
    plt.savefig('plots_review/figure3/bimodality_reduction_labels.png', bbox_inches='tight')
    plt.show()
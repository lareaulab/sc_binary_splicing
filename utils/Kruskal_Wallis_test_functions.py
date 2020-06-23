current_palette = sns.color_palette('muted')

from scipy.stats import friedmanchisquare
from scipy.stats import zscore
from tqdm import tqdm_notebook as tqdm
from scipy.special import logit

##########
def hyper_test(M, n, N, k):
    '''
    Calculates the hypergeometric test.
    
    Input:
      M: Population size (total exons tested)
      n: Successess in population (exons with p-value <= x)
      N: Sample size (exons selected)
      k: Successes in sample (selected exons with p-value <= x)
    Output:
      p_depleted: p-value of depletion from the hypergeometric test
      p_enriched: p-value of enrichment from the hypergeometric test
    '''
    hpd = hypergeom(M, n, N) # Hypergeometric distribution
    p_depleted = hpd.cdf(k)  # Get cummulative distribution up to k
    p_enriched = hpd.sf(k-1) # Survival function up to k-1 (1 - cdf(k-1))
    return p_depleted, p_enriched
    
    
##########
def linearize_psi(sample):
    '''
    Calculates the logit of PSI, but caps at 0.01 and 0.99
    
    Input:
      sample: a list of PSI values
    Output:
      linear_sample: logit of the PSI
    '''
    linear_sample = []
    for x in sample:
        if x >= 0.99:
            linear_sample.append(logit(0.99))
        elif x <= (0.01):
            linear_sample.append(logit(0.01))
        else:
            linear_sample.append(logit(x))
    return linear_sample
    
    
##########
def run_anova(samples, linearize = False):
    '''
    Runs the Kruskal-Wallis analysis of variance for data between 2 up to 6 clusters.
    This function is necessary because for some reason, the Python implementation of the
    Kruskal-Wallis test can take mutiple observations, but not as an array, thus I have 
    to resort to this function to test different number of clusters.
    
    Input:
      samples: list of list of PSI observations
    Output:
      returns Kruskal-Wallis test results over the list of lists
    '''
    
    if linearize:
        samples = [linearize_psi(x) for x in samples]
        samples = [x - np.mean(x) for x in samples]
        
    if len(samples) == 2:
        return kruskal(samples[0], samples[1])
    elif len(samples) == 3:
        return kruskal(samples[0], samples[1], samples[2])
    elif len(samples) == 4:
        return kruskal(samples[0], samples[1], samples[2], samples[3])
    elif len(samples) == 5:
        return kruskal(samples[0], samples[1], samples[2], samples[3], samples[4])
    elif len(samples) == 6:
        return kruskal(samples[0], samples[1], samples[2], samples[3], samples[4], samples[5])
    

def test_exon_bimodal_anova(PSI_tab, exon, pca_clust, clusters = 'AC', obs_min=0.5, linearize=False, n=3):
    '''
    Run Kruskal-Wallis test for one exon, to get significance in the differences 
    in PSI between multiple clusters. 
    
    Input
      PSI_tab: Matrix of PSI
      exon: name of the exon to test
      pca_clust: metadata dataframe with cluster information
      clusters: column of pca_clust with cluster information
      obs_min: minimum % of cells in cluster that have an observation to run the test (default: 50%)
      linearize: if calculate logit of the PSI (default: False)
    Output:
      anova_p: p-value of the KW test on the input exon
      10: vestigial. Downstream code expects two outputs
    '''
    
    obs_cells = PSI_tab.loc[exon].dropna().index # select cells that have a PSI observation
    cluster_psi = [] # list of lists of PSI observations per cluster
    
    for i in pca_clust[clusters].unique(): # iterate through each cluster
        clust_cells = pca_clust.loc[pca_clust[clusters] == i].index # select the cells in the cluster
        c_cells = [x for x in obs_cells if x in clust_cells] # cells in cluster with PSI observation of the target exon
        if len(c_cells)/len(clust_cells) >= obs_min: # Make sure that minimum % of cells in the cluster is met
            psi = list(PSI_tab.loc[exon, c_cells]) # list of PSI observations of the exon in cluster
            cluster_psi.append(psi)
    if len(cluster_psi) >= n: # run the test only if the exon is observed in at least n different clusters
        try:
            anova_p = run_anova(cluster_psi, linearize)[1] # p-value of the KW test
        except:
            anova_p = np.nan # if the test crashes for some reason
    else:
        anova_p = np.nan # not enough observations to run test
    return anova_p, len(cluster_psi), 10




def cluster_anova_test(PSI_tab, pca_clust, clusters, correction = 'fdr_bh', 
                          correct_multitest = True, obs_min = 0.5, linearize=False, n=3):
    '''
    Runs the Kruskal-Wallis test for a PSI matrix, and a given set of clusters.
    It wraps the test_exon_bimodal_anova function for all exons.
    
    Input
      PSI_tab: Matrix of PSI
      pca_clust: metadata dataframe with cluster information
      clusters: column of pca_clust with cluster information
      correction: vestigial; used to include an option to correct p-values
      correct_multitest: vestigial. p-values are not corrected anymore, because the significance of
                         individual exons is not the focus of this test.
      obs_min: minimum % of cells in cluster that have an observation to run the test 
               (default: 50%; for three clusters min)
      linearize: if calculate logit of the PSI (default: False)
    Output
      cluster_df: dataframe with p-values for each exon that meets the observed minimum
    '''
    cluster_array = []
    pvals = []
    exon_pass = []
    not_pass = 0
    
    for i in tqdm(range(len(PSI_tab.index))):
        
        exon = PSI_tab.index[i]
        anv_p, pos, neg = test_exon_bimodal_anova(PSI_tab, exon, pca_clust, clusters = 'AC', obs_min=obs_min, linearize=linearize, n=n)
        if not np.isnan(anv_p):
            cluster_array.append(pos/(neg+pos))
            pvals.append(anv_p)
            exon_pass.append(exon)
    

        else:
            not_pass += 1
        
    print('not pass: '+str(not_pass))
    print('tested exons: ' + str(len(pvals)))
    figsize(4, 4)
    plt.hist(pvals)
    plt.xlabel('p-value')
    plt.ylabel('occurrences')
    plt.show()
    
    cluster_df = pd.DataFrame()
    cluster_df['pvals'] = pvals
    cluster_df.index = exon_pass
    
    return cluster_df


def test_anova_filters(PSI_tab, mrna_counts, mrna_per_event, read_counts, coverage_tab, pca_clust, clusters = 'AC',
                     psi_min = 0.2, obs_min = 0.5, mrna_min = 10, mrna_read_min=0, read_min = 10, filter_obs = False, 
                    dset_name = '', correct_multitest = False, linearize=False, n=3):
    '''
    Wrapper function that manages the run of the Kruskal-Wallis test in the dataset, in addition to
    running basic filtering and exon selection functions. At the moment of writing this note, many parts 
    of the code are vestigial and will be removed.
    
    Input:
      PSI_tab: Matrix of PSI
      mrna_counts: matrix of mRNA molecules per gene
      mrna_per_event: mrna_counts with PSI_tab index; extended for multiple exons per gene
      read_counts: SJ counts used to estimate observations in PSI_tab
      coverage_tab: splice junction coverage rate
      pca_clust: metadata matrix with cluster information for cells
      clusters: column in pca_clust with clusters (default AC, but cell_type can also be used)
      psi_min: consider only exons with PSI average between [psi_min, 1-psi_min]
      obs_min: minimum % of cells in cluster that have an observation to run the test 
               (default: 50%; for three clusters min)
      mrna_min: minimum number of mRNAs in a quality observation (default 10)
      mrna_read_min: set an additional baseline minimum of reads for the mRNA filter (default 0)
      read_min: flat minimum of informative reads per event for the read filter (default 10)
      filter_obs: vestigial; should not be run as True
      dset_name: vestigial; this function used to make plots as well
      correct_multitest: vestigial
      linearize: whether if linearize the PSI before running KW test. We've found that it has
                 very little effect in the results.
    Output
      change_tab: table with KW p-values for exons with minimum observations
      mrna_selected: set of exons that pass the mRNA filter
      mrna_only_selected: vestigial; set of exons that pass the 10 mRNA in gene minimum, but not
                          necessarily the minimum reads expected for 10 mRNAs
      read_selected: set of exons that pass the flat-read minimum filter
    '''
    
    observed = PSI_tab.loc[PSI_tab.isna().mean(axis=1) <= (1-obs_min)].index
    int_genes, int_exons = spu.get_int_events(PSI_tab.loc[observed], mrna_counts, psi_min)

    aver = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, mrna_min, mrna_read_min, obs_min)
    
    aver_mrna_only = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, mrna_min, 0, obs_min, filter_cj = False)
    
    aver_reads = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, 0, read_min, obs_min)
    
    aver_all = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, mrna_min, mrna_read_min, 0)
    
    aver_all_mrna_only = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, mrna_min, 0, 0, filter_cj = False)
    
    aver_all_reads = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, 0, read_min, 0)

    joint_idx = [x for x in int_exons if x in aver_all[0].index]
    joint_idx = [x for x in joint_idx if x in aver_all_reads[0].index]
    joint_idx = [x for x in joint_idx if x in aver_all_mrna_only[0].index]
    print('intermediate exons: ' + str(len(joint_idx)))
    
    if filter_obs:
        change_tab = cluster_anova_test(aver_all[0].loc[joint_idx], pca_clust, clusters, obs_min=obs_min, linearize=linearize, n=n)
        change_tab_reads = cluster_anova_test(aver_all_reads[0].loc[joint_idx], pca_clust, clusters, obs_min=obs_min, linearize=linearize, n=n)
    else:
        change_tab = cluster_anova_test(PSI_tab.loc[joint_idx], pca_clust, clusters, obs_min=obs_min, linearize=linearize, n=n)
        
    if correct_multitest:
        adj_pvals = multipletests(list(change_tab.pvals), method='fdr_bh')[1]
        change_tab['pvals'] = list(adj_pvals)

    joint_idx = change_tab.index
    mrna_selected = aver[0].index
    mrna_only_selected = aver_mrna_only[0].index
    read_selected = aver_reads[0].index
    mrna_selected = list([x for x in joint_idx if x in mrna_selected])
    mrna_only_selected = list([x for x in joint_idx if x in mrna_only_selected])
    read_selected = list([x for x in joint_idx if x in read_selected])
    
    return change_tab, mrna_selected, mrna_only_selected, read_selected


###########
###########
# The following functions run a similar test, but using Spearman correlation with pseudotime
# It is currently not used for anything related to the Buen Abad Najar et al., 2019 pre-print.

def test_exon_bimodal_spearman(PSI_tab, exon, obs_min, pca_clust, clusters = 'AC', one_test = False, linearize=False):
    
    obs_cells = PSI_tab.loc[exon].dropna().index
    
    if one_test:
        
        psi = list(PSI_tab.loc[exon, obs_cells])
        
        if linearize:
            psi = linearize_psi(psi)
        
        pseudotime =  list(pca_clust.loc[obs_cells, 'pseudotime'])
        spearman_p = spearmanr(pseudotime, psi)[1]
        
    else:
        
        pvals = []
        
        for i in range(len(pca_clust[clusters].unique()) - 1):
            
            c1 = i
            c2 = i+1
            c1_ = pca_clust.loc[pca_clust[clusters] == c1].index
            c2_ = pca_clust.loc[pca_clust[clusters] == c2].index
            c1_cells = [x for x in obs_cells if x in c1_]
            c2_cells = [x for x in obs_cells if x in c2_]
            
            cluster_cells = list(c1_cells) + list(c2_cells)
            psi = list(PSI_tab.loc[exon, cluster_cells])
            
            if linearize:
                psi = linearize_psi(psi)
            
            pseudotime =  list(pca_clust.loc[cluster_cells, 'pseudotime'])
            
            if ((len(c1_cells)/len(c1_)) >= obs_min) and ((len(c1_cells)/len(c1_)) >= obs_min):
                pvals.append(spearmanr(pseudotime, psi)[1])
                
        if len(pvals) >= 2:
            spearman_p = combine_pvalues(pvals)[1]
        else:
            spearman_p = np.nan
        
    return spearman_p, np.mean(pseudotime), np.mean(psi)


def cluster_spearman_test(PSI_tab, pca_clust, clusters, obs_min, correction = 'fdr_bh', 
                          correct_multitest = True, one_test = False, linearize=False):
    cluster_array = []
    pvals = []
    exon_pass = []
    not_pass = 0
    for i in tqdm(range(len(PSI_tab.index))):

        exon = PSI_tab.index[i]
        chi_p, pos, neg = test_exon_bimodal_spearman(PSI_tab, exon, obs_min, pca_clust, clusters = 'AC', 
                                                     one_test = one_test, linearize=linearize)
        if not np.isnan(chi_p):
            cluster_array.append(pos/(neg+pos))
            pvals.append(chi_p)
            exon_pass.append(exon)
            

        else:
            not_pass += 1
        
    print('not pass')
    print(not_pass)
    figsize(4, 4)
    plt.hist(pvals)
    plt.show()
    
    cluster_df = pd.DataFrame()
    cluster_df['pvals'] = pvals
    cluster_df.index = exon_pass
    
    return cluster_df


def test_spearman_filters(PSI_tab, mrna_counts, mrna_per_event, read_counts, coverage_tab, pca_clust, clusters = 'AC',
                     psi_min = 0.2, obs_min = 0.5, mrna_min = 10, mrna_read_min=10, read_min = 10, filter_obs = False, 
                    dset_name = '', one_test = False, correct_multitest=False, linearize=False):
    
    observed = PSI_tab.loc[PSI_tab.isna().mean(axis=1) <= (1-obs_min)].index
    int_genes, int_exons = spu.get_int_events(PSI_tab.loc[observed], mrna_counts, psi_min)


    aver = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, mrna_min, mrna_read_min, obs_min)
    
    aver_mrna_only = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, mrna_min, 0, obs_min, filter_cj=False)
    
    aver_reads = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, 0, read_min, obs_min)
    
    aver_all = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, mrna_min, mrna_read_min, 0)
    
    aver_all_mrna_only = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, mrna_min, 0, 0, filter_cj=False)
    
    aver_all_reads = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, 0, read_min, 0)

    joint_idx = [x for x in int_exons if x in aver_all[0].index]
    joint_idx = [x for x in joint_idx if x in aver_all_reads[0].index]
    joint_idx = [x for x in joint_idx if x in aver_all_mrna_only[0].index]
    print('intermediate exons: ' + str(len(joint_idx)))
    
    if filter_obs:
        change_tab, chi_p = cluster_spearman_test(aver_all[0].loc[joint_idx], pca_clust, clusters, 
                                                  one_test = one_test, linearize=linearize)
        change_tab_reads, chi_p_reads = cluster_spearman_test(aver_all_reads[0].loc[joint_idx], pca_clust, 
                                                              clusters, one_test = one_test, linearize=linearize)
        change_tab_reads['pvals'] = chi_p_reads
    else:
        change_tab = cluster_spearman_test(PSI_tab.loc[joint_idx], pca_clust, clusters, obs_min, one_test = one_test, linearize=linearize)
        
    
    if correct_multitest:
        adj_pvals = multipletests(list(change_tab.pvals), method='fdr_bh')[1]
        change_tab['pvals'] = list(adj_pvals)

    joint_idx = change_tab.index
    mrna_selected = aver[0].index
    mrna_only_selected = aver_mrna_only[0].index
    read_selected = aver_reads[0].index
    mrna_selected = list([x for x in joint_idx if x in mrna_selected])
    mrna_only_selected = list([x for x in joint_idx if x in mrna_only_selected])
    read_selected = list([x for x in joint_idx if x in read_selected])
    
    return change_tab, mrna_selected, mrna_only_selected, read_selected

########
########

def get_statistics(pvals, selected_exons, x, beta):
    '''
    Get several summary calculations to evaluate the performance of an exon selection
    approach. These are inspired in machine learning classifier methods. 
    
    It basically treats the exon selection approach as a classifier that predicts
    exons with significant splicing changes. The true instances are based in a minimum
    p-value cutoff x. We set:
    
    For a given p-value cutoff x: 
    — The set of “True instances” X1 is all exons with p-value <= x
    — The set of “False instances” X2 is all exons with p-value > x
    — The set of “Positive calls”  X3 is all exons selected by the filter
    — The set of “Negative calls”  X4 is all exons rejected by the filter

    TP = X1 AND X3
    TN = X2 AND X4
    FP = X2 AND X3
    FN = X1 AND X4
    '''

    X1 = set(pvals.loc[np.array(pvals.pvals) <= x].index)
    X2 = set(pvals.loc[np.array(pvals.pvals) > x].index)
    X3 = set(selected_exons)
    X4 = set([x for x in pvals.index if x not in selected_exons])

    TP = X1 & X3
    TN = X2 & X4
    FP = X2 & X3
    FN = X1 & X4

    try:
        recall = len(TP)/len(X1) # Also called sensitivity, or TPR
    except:
        recall = np.nan
    precision = len(TP)/len(X3)
    specificity = len(TN)/len(X2)
    try:
        F1_score = (1+(beta**2))*(precision*recall)/(((beta**2)*precision)+recall)
    except:
        F1_score = np.nan
    hyper_recall_specificity = (1+(beta**2))*(specificity*recall)/(((beta**2)*specificity)+recall)
    
    try:
        DOR = (recall*specificity)/((1-recall)*(1-specificity))
    except:
        DOR = np.nan
        
    ACC = (len(TP) + len(TN)) / (len(X1)+len(X2))
    BA = (recall*specificity)/2
    
    return recall, precision, specificity, F1_score, hyper_recall_specificity, DOR, ACC, BA


def summary_curves(pvals, selected_exons, p_limits, beta):
    recall_list = []
    precision_list = []
    specificity_list = []
    f1_score_list = []
    hrs_list = []
    DOR_list = []
    ACC_list = []
    BA_list = []
    
    for x in p_limits:
        
        recall, precision, specificity, f1_score, hrs, DOR , ACC, BA = get_statistics(pvals, selected_exons, x, beta)
        recall_list.append(recall)
        precision_list.append(precision)
        specificity_list.append(specificity)
        f1_score_list.append(f1_score)
        hrs_list.append(hrs)
        DOR_list.append(DOR)
        ACC_list.append(ACC)
        BA_list.append(BA)
        
    return recall_list, precision_list, specificity_list, f1_score_list, hrs_list, DOR_list, ACC_list, BA_list
    


def summary_plots(data_out, dset_name, 
                  p_lim=1, p_low = 0.0001, p_steps = 0.0001, linear_p = True, beta=1, plot_line=True):

    pvals, selected_mrna, selected_mrna_only, selected_read = data_out
    
    if linear_p:
        p_limits = np.arange(p_low, p_lim + p_steps, p_steps)
        
    else:
        p_limits = 10**(-np.arange(-np.log10(p_lim), -np.log10(p_low)+p_steps, p_steps))
    
    recall_mrna, precision_mrna, specificity_mrna, f1_score_mrna, hrs_mrna, DOR_mrna, ACC_mrna, BA_mrna = summary_curves(pvals, selected_mrna, p_limits, beta)
    recall_mrna_only, precision_mrna_only, specificity_mrna_only, f1_score_mrna_only, hrs_mrna_only, DOR_mrna_only, ACC_mrna_only, BA_mrna_only = summary_curves(pvals, selected_mrna_only, p_limits, beta)
    recall_read, precision_read, specificity_read, f1_score_read, hrs_read, DOR_read, ACC_read, BA_read = summary_curves(pvals, selected_read, p_limits, beta)
    
    
    figsize(6,4)
    fig = plt.figure()
    ax  = plt.subplot(1,1,1)
    
    ax.plot(-np.log10(p_limits), recall_read, c='darkred', label='read filter', linewidth = 3)
    ax.plot(-np.log10(p_limits), recall_mrna, c='navy', label='mRNA filter', linewidth = 3)
    
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(labelsize=28, length=5)
    plt.title('Recall curves, ' + dset_name, fontsize=28)
    plt.xlabel('-log10 p-value', fontsize=28)
    plt.ylabel('Recall', fontsize=28)
    plt.xticks([1, 2, 3, 4, 5], ['1', '2', '3', '4', '5'])
    plt.show()
    
    
    fig = plt.figure()
    ax  = plt.subplot(1,1,1)
    
    ax.plot(-np.log10(p_limits), precision_read, c='darkred', label='read filter', linewidth = 3)
    ax.plot(-np.log10(p_limits), precision_mrna, c='navy', label='mRNA filter', linewidth=3)
    
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(labelsize=28, length=5)
    plt.title('Precision curves, ' + dset_name, fontsize=28)
    plt.xlabel('-log10 p-value', fontsize=28)
    plt.ylabel('Precision', fontsize=28)
    plt.xticks([1, 2, 3, 4, 5], ['1', '2', '3', '4', '5'])
    
    plt.show()
    
    fig = plt.figure()
    ax  = plt.subplot(1,1,1)
    
    ax.plot(-np.log10(p_limits), specificity_read, c='darkred', label='read filter', linewidth=3)
    ax.plot(-np.log10(p_limits), specificity_mrna, c='navy', label='mRNA filter', linewidth=3)
    
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(labelsize=28, length=5)
    plt.title('Specificity curves, ' + dset_name, fontsize=28)
    plt.xlabel('-log10 p-value', fontsize=28)
    plt.ylabel('Specificity', fontsize=28)
    plt.xticks([1, 2, 3, 4, 5], ['1', '2', '3', '4', '5'])
    
    plt.show()
    
    fig = plt.figure()
    ax  = plt.subplot(1,1,1)
    
    ax.plot(-np.log10(p_limits), f1_score_read, c='darkred', label='read filter', linewidth=3)
    ax.plot(-np.log10(p_limits), f1_score_mrna, c='navy', label='mRNA filter', linewidth=3)

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(labelsize=28, length=5)
    plt.title('F'+str(round(beta, 2))+' curves, ' + dset_name, fontsize=28)
    plt.xlabel('-log10 p-value', fontsize=28)
    plt.ylabel('F'+str(round(beta, 2))+' score', fontsize=28)
    plt.xticks([1, 2, 3, 4, 5], ['1', '2', '3', '4', '5'])
    
    plt.show()
    
    
    fig = plt.figure()
    ax  = plt.subplot(1,1,1)
    
    ax.plot(-np.log10(p_limits), hrs_read, c='darkred', label='read filter', linewidth=3)
    ax.plot(-np.log10(p_limits), hrs_mrna, c='navy', label='mRNA filter', linewidth=3)
    
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(labelsize=28, length=5)
    plt.title('Sensitivity/specificity, ' + dset_name, fontsize=28)
    plt.xlabel('-log10 p-value', fontsize=28)
    plt.ylabel('Harmonic mean', fontsize=28)
    plt.xticks([1, 2, 3, 4, 5], ['1', '2', '3', '4', '5'])

    plt.show()
    
    
    fig = plt.figure()
    ax  = plt.subplot(1,1,1)
    
    ax.plot(-np.log10(p_limits), DOR_read, c='darkred', label='read filter', linewidth=3)
    ax.plot(-np.log10(p_limits), DOR_mrna, c='navy', label='mRNA filter', linewidth=3)
    
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(labelsize=28, length=5)
    plt.title('Diagnostic odds ratio, ' + dset_name, fontsize=28)
    plt.xlabel('-log10 p-value', fontsize=28)
    plt.ylabel('DOR', fontsize=28)
    plt.xticks([1, 2, 3, 4, 5], ['1', '2', '3', '4', '5'])
    
    plt.show()
    

    fig = plt.figure()
    ax  = plt.subplot(1,1,1)
    
    ax.plot(-np.log10(p_limits), ACC_read, c='darkred', label='read filter', linewidth=3)
    ax.plot(-np.log10(p_limits), ACC_mrna, c='navy', label='mRNA filter', linewidth=3)
    
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(labelsize=28, length=5)
    plt.title('Accuracy, ' + dset_name, fontsize=28)
    plt.xlabel('-log10 p-value', fontsize=28)
    plt.ylabel('ACC', fontsize=28)
    plt.xticks([1, 2, 3, 4, 5], ['1', '2', '3', '4', '5'])

    plt.show()
    
    
    fig = plt.figure()
    ax  = plt.subplot(1,1,1)
    
    ax.plot(-np.log10(p_limits), BA_read, c='darkred', label='read filter', linewidth=3)
    ax.plot(-np.log10(p_limits), BA_mrna, c='navy', label='mRNA filter', linewidth=3)
    
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(labelsize=28, length=5)
    plt.title('Balanced Accuracy, ' + dset_name, fontsize=28)
    plt.xlabel('-log10 p-value', fontsize=28)
    plt.ylabel('BA', fontsize=28)
    plt.xticks([1, 2, 3, 4, 5], ['1', '2', '3', '4', '5'])
    
    plt.show()
    
    
    
    fig = plt.figure()
    ax  = plt.subplot(1,1,1)
    
    linea = []
    for x in p_limits:
        linea.append(np.sum(pvals <= x))
    
    ax.plot(-np.log10(p_limits), linea, c='navy', label='total exons', linewidth=3)

    
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(labelsize=28, length=5)
    plt.title('Significant exons, ' + dset_name, fontsize=28)
    plt.xlabel('-log10 p-value', fontsize=28)
    plt.ylabel('Total significant exons', fontsize=28)
    plt.xticks([1, 2, 3, 4, 5], ['1', '2', '3', '4', '5'])

    plt.show()

    
    figsize(4,3)
    fig = plt.figure()
    ax  = plt.subplot(1,1,1)

    ax.plot([], c='darkred', label='read filter', linewidth=3)
    ax.plot([], c='forestgreen', label='mRNA only filter', linewidth=3)
    ax.plot([], c='navy', label='mRNA filter', linewidth=3)

    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)
    ax.xaxis.set_ticks_position('none') 
    ax.yaxis.set_ticks_position('none')
    plt.legend(fontsize=28, frameon=False)

    plt.show()

    
    
    
def plot_pvalues(ax, test_data):
    pvals, selected_mrna, selected_mrna_only, selected_read = test_data
    current_palette = sns.color_palette('dark')
    pplot_1 = probplot(pvals.pvals, dist='uniform')[0]
    pplot_2 = probplot(pvals.loc[selected_read, 'pvals'], dist='uniform')[0]
    pplot_4 = probplot(pvals.loc[selected_mrna, 'pvals'], dist='uniform')[0]
    ax.plot(pplot_1[0], pplot_1[1], alpha = 0.75, c=current_palette[4], label='all exons', linewidth=2)
    ax.plot(pplot_2[0], pplot_2[1], alpha = 0.75, c='darkred', label='read filter', linewidth=2)
    ax.plot(pplot_4[0], pplot_4[1], alpha = 0.75, c='navy', label='combined filter', linewidth=2)
    ax.plot([0, 1], [0.05, 0.05], 'r--', linewidth=1)
    ax.set_xlabel('Quantiles', fontsize=28)
    ax.tick_params(labelsize=28, length=5)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels(['0', '0.5', '1'])
    ax.tick_params(labelsize=28, length=5)
    
    
    
def get_fold_enrichment(test_data, pval):
    
    pvals, selected_mrna, selected_mrna_only, selected_read = test_data
    M = len(pvals)
    n = np.sum((pvals <= pval).pvals)
    
    N_mrna = len(selected_mrna)
    k_mrna = np.sum((pvals.loc[selected_mrna] <= pval).pvals)
    
    N_mrna_only = len(selected_mrna_only)
    k_mrna_only = np.sum((pvals.loc[selected_mrna_only] <= pval).pvals)
    
    N_read = len(selected_read)
    k_read = np.sum((pvals.loc[selected_read] <= pval).pvals)
    
    mrna_enrichment = (k_mrna*M)/(n*N_mrna)
    mrna_only_enrichment = (k_mrna_only*M)/(n*N_mrna_only)
    read_enrichment = (k_read*M)/(n*N_read)
    
    return read_enrichment, mrna_only_enrichment, mrna_enrichment, M, n, N_read, k_read, N_mrna_only, k_mrna_only, N_mrna, k_mrna 
    
def get_enrichment(test_data):
    p_limits = 10**(-np.arange(-np.log10(0.1), -np.log10(0.00001)+0.05, 0.05))
    read_enrichment = []
    mrna_only_enrichment = []
    mrna_enrichment = []
    
    read_pvals = []
    mrna_only_pvals = []
    mrna_pvals = []
    
    for x in p_limits:
        enrich = get_fold_enrichment(test_data, x)
        read_enrichment.append(enrich[0])
        mrna_only_enrichment.append(enrich[1])
        mrna_enrichment.append(enrich[2])
        
        read_pvals.append(hyper_test(enrich[3], enrich[4], enrich[5], enrich[6])[1])
        mrna_only_pvals.append(hyper_test(enrich[3], enrich[4], enrich[7], enrich[8])[1])
        mrna_pvals.append(hyper_test(enrich[3], enrich[4], enrich[9], enrich[10])[1])
        
    read_padj = multipletests(read_pvals, method='fdr_bh')[1]
    mrna_only_padj = multipletests(mrna_only_pvals, method='fdr_bh')[1]
    mrna_padj = multipletests(mrna_pvals, method='fdr_bh')[1]

#     corrected = multipletests(read_pvals + mrna_pvals, method='fdr_bh')[1]
#     read_padj = corrected[:len(p_limits)]
# #     mrna_only_padj = corrected[len(p_limits):2*len(p_limits)]
#     mrna_only_padj = corrected[len(p_limits):]
#     mrna_padj = corrected[len(p_limits):]
        
    return read_enrichment, mrna_only_enrichment, mrna_enrichment, read_padj, mrna_only_padj, mrna_padj


def plot_filter_lines(ax, read, mrna_only, mrna):
    current_palette = sns.color_palette('dark')
    ax.plot((np.arange(-np.log10(0.1), -np.log10(0.00001)+0.05, 0.05)), read, c='darkred', 
            label='read filter', linewidth=2)
    ax.plot((np.arange(-np.log10(0.1), -np.log10(0.00001)+0.05, 0.05)), mrna, c='navy', 
            label='combined filter', linewidth=2)
    ax.set_xlabel('-log10 p-value', fontsize=28)
    ax.tick_params(labelsize=28, length=5)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set_xticks([1, 2, 3, 4, 5])
    ax.set_xticklabels(['1', '2', '3', '4', '5'])
    ax.tick_params(labelsize=28, length=5)
    
    
    
    
def plot_curves(ax, test_data, which_curve = 0, p_max = 0.1, p_min = 0.00001, steps = 0.05):
    
    names = ['Recall', 'Precision', 'Specificity', 'F1 score', 'Recall/Specificity', 
             'LR+', 'LR-', 'Diagnostic Odds Ratio', 'Accuracy', 'Balanced Accuracy']
    
    name = names[which_curve]    
    pvals, selected_mrna, selected_read = test_data
    p_limits = 10**(-np.arange(-np.log10(p_max), -np.log10(p_min)+steps, steps))
    summary_mrna = summary_curves(pvals, selected_mrna, p_limits, 1)[which_curve]
    summary_read = summary_curves(pvals, selected_read, p_limits, 1)[which_curve]
    ax.plot(-np.log10(p_limits), summary_read, c='darkred', label='read filter', linewidth=2)
    ax.plot(-np.log10(p_limits), summary_mrna, c='navy', label='mRNA filter', linewidth=2)
    ax.set_xlabel('-log10 p-value', fontsize=28)
    ax.tick_params(labelsize=28, length=5)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    
    ax.tick_params(labelsize=28, length=5)

    
def plot_dset_curves(chen_data_out, trapnell_data_out, song_data_out, 
                     which_curve = 0, p_max = 0.1, p_min = 0.00001, steps = 0.05, ylims=(0, 1.05)):
    
    names = ['Recall', 'Precision', 'Specificity', 'F1 score', 'Recall/Specificity', 
             'LR+', 'LR-', 'Diagnostic Odds Ratio', 'Accuracy', 'Balanced Accuracy']
    
    figsize(18, 5)
    fig = plt.figure()
    gs = GridSpec(1,3)
    gs.update(wspace=0.05, hspace=1.05)
    ax_1 = fig.add_subplot(gs[0,0])
    ax_2 = fig.add_subplot(gs[0,1])
    ax_3 = fig.add_subplot(gs[0,2])
    ax_1.set_ylim(ylims)
    ax_2.set_ylim(ylims)
    ax_3.set_ylim(ylims)
    plot_curves(ax_1, chen_data_out, which_curve, p_max, p_min, steps)
    plot_curves(ax_2, trapnell_data_out, which_curve, p_max, p_min, steps)
    plot_curves(ax_3, song_data_out, which_curve, p_max, p_min, steps)
    xticks = range(int(-np.log10(p_max)), int(-np.log10(p_min))+1)
    ax_1.set_xticks(xticks)
    ax_1.set_xticklabels([str(i) for i in xticks])
    ax_2.set_xticks(xticks)
    ax_2.set_xticklabels([str(i) for i in xticks])
    ax_3.set_xticks(xticks)
    ax_3.set_xticklabels([str(i) for i in xticks])
    ax_1.set_ylabel(names[which_curve], fontsize = 28)
    plt.setp(ax_2.get_yticklabels(), visible=False)
    plt.setp(ax_3.get_yticklabels(), visible=False)
    ax_2.yaxis.set_ticks_position('none')
    ax_3.yaxis.set_ticks_position('none')    
    plot_name = '_'.join(names[which_curve].split()).lower()

    plt.savefig('plots_review/figure3/anova/' + plot_name + '.svg', bbox_inches='tight', transparent=True)
    plt.savefig('plots_review/figure3/anova/' + plot_name + '.pdf', bbox_inches='tight', transparent=True)
    plt.savefig('plots_review/figure3/anova/' + plot_name + '.png', dpi=300, bbox_inches='tight', transparent=True)

    plt.show()

    
    
def run_hypergeom_test(test_data):
    pvals, selected_mrna, selected_mrna_only, selected_read = test_data
    M = len(pvals)
    n = np.sum((pvals <= 0.05).pvals)
    
    N_mrna = len(selected_mrna)
    k_mrna = np.sum((pvals.loc[selected_mrna] <= 0.05).pvals)
    
    N_mrna_only = len(selected_mrna_only)
    k_mrna_only = np.sum((pvals.loc[selected_mrna_only] <= 0.05).pvals)
    
    N_read = len(selected_read)
    k_read = np.sum((pvals.loc[selected_read] <= 0.05).pvals)
    
    hyper_pvals = []
    hyper_pvals.append(hyper_test(M, n, N_read, k_read)[1])
    hyper_pvals.append(hyper_test(M, n, N_mrna_only, k_mrna_only)[1])
    hyper_pvals.append(hyper_test(M, n, N_mrna, k_mrna)[1])
    
    
    return M, n, [N_read, N_mrna_only, N_mrna], [k_read, k_mrna_only, k_mrna], hyper_pvals
def hyper_test(M, n, N, k):
    
    hpd = hypergeom(M, n, N)
    p_depleted = hpd.cdf(k)
    p_enriched = hpd.sf(k-1)
    
    return p_depleted, p_enriched
    
    
def run_anova(samples):
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
    

def test_exon_bimodal_anova(PSI_tab, exon, pca_clust, clusters = 'AC', obs_min=0.5):
    
    obs_cells = PSI_tab.loc[exon].dropna().index
    
    cluster_psi = []
    
    for i in pca_clust[clusters].unique():
#         t1 = time.time()
        clust_cells = pca_clust.loc[pca_clust[clusters] == i].index
        c_cells = [x for x in obs_cells if x in clust_cells]
#         print('time1')
#         print(time.time()-t1)
#         if len(c_cells) >= 10: #########################
        if len(c_cells)/len(clust_cells) >= obs_min:
            psi = list(PSI_tab.loc[exon, c_cells])
            cluster_psi.append(psi)
    if len(cluster_psi) >= 3:
        anova_p = run_anova(cluster_psi)[1]
    else:
        anova_p = np.nan
        
        
    return anova_p, len(cluster_psi), 10

from tqdm import tqdm_notebook as tqdm

def cluster_anova_test(PSI_tab, pca_clust, clusters, correction = 'fdr_bh', 
                          correct_multitest = True, obs_min = 0.5):
    cluster_array = []
    pvals = []
    exon_pass = []
    not_pass = 0
#     completed = 0
    
    for i in tqdm(range(len(PSI_tab.index))):
        
#       for exon in PSI_tab.index:
#       t = time.time()
        exon = PSI_tab.index[i]
        anv_p, pos, neg = test_exon_bimodal_anova(PSI_tab, exon, pca_clust, clusters = 'AC', obs_min=obs_min)
        if not np.isnan(anv_p):
            cluster_array.append(pos/(neg+pos))
            pvals.append(anv_p)
            exon_pass.append(exon)
    

        else:
            not_pass += 1
#         completed += 1
#         if completed % 100 == 0:
#             print (completed)
        
    print('not pass: '+str(not_pass))
    print('tested exons: ' + str(len(pvals)))
#     print(not_pass)
#     figsize(4, 4)
#     plt.hist(pvals)
#     plt.show()
#     if correct_multitest:
#         pvals_adj = multipletests(pvals, method='fdr_bh')[1]
#     else:
#         pvals_adj = pvals
    cluster_df = pd.DataFrame()
    cluster_df['pvals'] = pvals
#     cluster_df['array'] = cluster_array
    cluster_df.index = exon_pass
    
    return cluster_df


def test_anova_filters(PSI_tab, mrna_counts, mrna_per_event, read_counts, coverage_tab, pca_clust, clusters = 'AC',
                     psi_min = 0.2, obs_min = 0.5, mrna_min = 10, read_min = 10, filter_obs = False, 
                    dset_name = '', correct_multitest = False):
    
    observed = PSI_tab.loc[PSI_tab.isna().mean(axis=1) <= (1-obs_min)].index
    int_genes, int_exons = spu.get_int_events(PSI_tab.loc[observed], mrna_counts, psi_min)


    aver = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, mrna_min, 0, obs_min)
    
    aver_mrna_only = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, mrna_min, 0, obs_min, filter_cj = False)
    
    aver_reads = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, 0, read_min, obs_min)
    
    aver_all = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, mrna_min, 0, 0)
    
    aver_all_mrna_only = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, mrna_min, 0, 0, filter_cj = False)
    
    aver_all_reads = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, 0, read_min, 0)

    joint_idx = [x for x in int_exons if x in aver_all[0].index]
    joint_idx = [x for x in joint_idx if x in aver_all_reads[0].index]
    joint_idx = [x for x in joint_idx if x in aver_all_mrna_only[0].index]
    print('intermediate exons: ' + str(len(joint_idx)))
    
    if filter_obs:
        change_tab = cluster_anova_test(aver_all[0].loc[joint_idx], pca_clust, clusters, obs_min=obs_min)
        change_tab_reads = cluster_anova_test(aver_all_reads[0].loc[joint_idx], pca_clust, clusters, obs_min=obs_min)
#         change_tab_reads['pvals'] = chi_p_reads
    else:
        change_tab = cluster_anova_test(PSI_tab.loc[joint_idx], pca_clust, clusters, obs_min=obs_min)
        
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

def test_exon_bimodal_chi(PSI_tab, exon, clusters, psi_lim = 0.25):
    
    psi_mean = PSI_tab.loc[exon].mean()
    chi_pos = []
    chi_exp = []
    chi_neg = []
    for clust in clusters.unique():
        cluster_cells = clusters.loc[clusters == clust].index
        
        pos_cases = (PSI_tab.loc[exon, cluster_cells] >= (1-psi_lim)).sum()
        neg_cases = (PSI_tab.loc[exon, cluster_cells] <= psi_lim).sum()
        
        expected_pos = (pos_cases + neg_cases) * psi_mean
        
        
        chi_pos.append(pos_cases)
        chi_exp.append(expected_pos)
        chi_neg.append(neg_cases)
        
    chip = chisquare(chi_pos, chi_exp)[1]#
        
    return chip, np.array(chi_pos), np.array(chi_neg), np.array(chi_exp)#


def cluster_chi_test(PSI_tab, clusters, psi_lim, correction = 'fdr_bh', correct_multitest = True, 
                     print_extreme=False, extreme = 0.9):
    cluster_array = []
    pvals = []
    exon_pass = []
    not_pass = 0
    for i in tqdm(range(len(PSI_tab.index))):

        exon = PSI_tab.index[i]
        chi_p, pos, neg, exp = test_exon_bimodal_chi(PSI_tab, exon, clusters, psi_lim = psi_lim)
        if not np.isnan(chi_p):
            cluster_array.append(pos/(neg+pos))
            pvals.append(chi_p)
            exon_pass.append(exon)
            
            if print_extreme:
                if chi_p >= extreme:
                    print(pos)
                    print(neg)
                    print(exp)
        else:
            not_pass += 1
        
    print('not pass: '+str(not_pass))
    print('tested exons: ' + str(len(pvals)))
    
    cluster_df = pd.DataFrame()
    cluster_df['pvals'] = pvals
    cluster_df.index = exon_pass
    
    return cluster_df#, pvals


def test_chi_filters(PSI_tab, mrna_counts, mrna_per_event, read_counts, coverage_tab, clusters, 
                     psi_min = 0.2, psi_bin = 0.25, obs_min = 0.5, mrna_min = 10, read_min = 10, filter_obs = False, 
                    dset_name = '', correct_multitest = False):
    
    observed = PSI_tab.loc[PSI_tab.isna().mean(axis=1) <= (1-obs_min)].index
    int_genes, int_exons = spu.get_int_events(PSI_tab.loc[observed], mrna_counts, psi_min)


    aver = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, mrna_min, 0, obs_min)
    
    aver_mrna_only = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, mrna_min, 0, obs_min, filter_cj = False)
    
    aver_reads = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, 0, read_min, obs_min)
    
    aver_all = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, mrna_min, 0, 0)
    
    aver_all_mrna_only = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, mrna_min, 0, 0, filter_cj = False)
    
    aver_all_reads = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, 0, read_min, 0)

    joint_idx = [x for x in int_exons if x in aver_all[0].index]
    joint_idx = [x for x in joint_idx if x in aver_all_reads[0].index]
    joint_idx = [x for x in joint_idx if x in aver_all_mrna_only[0].index]
    print('intermediate exons: ' + str(len(joint_idx)))
    
    if filter_obs:
        change_tab, chi_p = cluster_chi_test(aver_all[0].loc[joint_idx], clusters, psi_bin)
        change_tab_reads, chi_p_reads = cluster_chi_test(aver_all_reads[0].loc[joint_idx], clusters, psi_bin)
        change_tab_reads['pvals'] = chi_p_reads
    else:
        change_tab = cluster_chi_test(PSI_tab.loc[joint_idx], clusters, psi_bin)
        

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

def test_exon_bimodal_hypergeom(PSI_tab, exon, clusters, psi_lim = 0.25, only_enrichment = True):
    
    psi_mean = PSI_tab.loc[exon].mean()
    
    run_test = True

    M = (~PSI_tab.loc[exon].isna()).sum()
    
    n_pos = (PSI_tab.loc[exon] >= (1-psi_lim)).sum()
    n_neg = (PSI_tab.loc[exon] <= psi_lim).sum()
    
    if (n_neg <= 10) or (n_pos <= 10):
#         print('reject, not enough pos or neg in general')
        run_test = False
    
    pvals = []
    hyp_pos = []
    hyp_neg = []
    
    ### Next: enrichment of only 1, or only 0; or both independently; related to all observations
    
    # Also, test just enrichment, vs enrichment + depletion
    
#     run_test = True
    
    
    for clust in clusters.unique():
        cluster_cells = clusters.loc[clusters == clust].index
        
        pos_cases = (PSI_tab.loc[exon, cluster_cells] >= (1-psi_lim)).sum()
        neg_cases = (PSI_tab.loc[exon, cluster_cells] <= psi_lim).sum()
        
        N = (~PSI_tab.loc[exon, cluster_cells].isna()).sum()
        
        if (pos_cases + neg_cases) < 1:
#             print('reject, small cluster')
            run_test = False
            pvals.append(np.nan)
        else:
            if only_enrichment:
                pvals.append(hyper_test(M, n_pos, N, pos_cases)[1])
                pvals.append(hyper_test(M, n_neg, N, neg_cases)[1])
            else:
                pvals.extend(list(hyper_test(M, n_pos, N, pos_cases)))
                pvals.extend(list(hyper_test(M, n_neg, N, neg_cases)))
            
        hyp_pos.append(pos_cases)
        hyp_neg.append(neg_cases)
        
    if run_test:
        hyper_p = combine_pvalues(pvals)[1]
    else:
        hyper_p = np.nan
        
#     if (hyper_p >= 0.99) or (hyper_p <= 0.01):
#         print(hyper_p)
#         print(M)
#         print(hyp_pos)
#         print(hyp_neg)
        
    return hyper_p, np.array(hyp_pos), np.array(hyp_neg)


def cluster_hypergeom_test(PSI_tab, clusters, psi_lim, correction = 'fdr_bh', correct_multitest = True):
    cluster_array = []
    pvals = []
    exon_pass = []
    not_pass = 0
    for i in tqdm(range(len(PSI_tab.index))):

        exon = PSI_tab.index[i]
        chi_p, pos, neg = test_exon_bimodal_hypergeom(PSI_tab, exon, clusters, psi_lim = psi_lim)
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
#     if correct_multitest:
#         pvals_adj = multipletests(pvals, method=correction)[1]
#     else:
#         pvals_adj = pvals
    cluster_df = pd.DataFrame()
    cluster_df['pvals'] = pvals
    cluster_df.index = exon_pass
    
    return cluster_df#, pvals


def test_hypergeom_filters(PSI_tab, mrna_counts, mrna_per_event, read_counts, coverage_tab, clusters, 
                     psi_min = 0.2, psi_bin = 0.25, obs_min = 0.5, mrna_min = 10, read_min = 10, filter_obs = False, 
                    dset_name = '', correct_multitest=False):
    
    observed = PSI_tab.loc[PSI_tab.isna().mean(axis=1) <= (1-obs_min)].index
    int_genes, int_exons = spu.get_int_events(PSI_tab.loc[observed], mrna_counts, psi_min)


    aver = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, mrna_min, 0, obs_min)
    
    aver_mrna_only = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, mrna_min, 0, obs_min, filter_cj = False)
    
    aver_reads = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, 0, read_min, obs_min)
    
    aver_all = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, mrna_min, 0, 0)
    
    aver_all_mrna_only = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, mrna_min, 0, 0, filter_cj = False)
    
    aver_all_reads = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, 0, read_min, 0)

    joint_idx = [x for x in int_exons if x in aver_all[0].index]
    joint_idx = [x for x in joint_idx if x in aver_all_reads[0].index]
    joint_idx = [x for x in joint_idx if x in aver_all_mrna_only[0].index]
    print('intermediate exons: ' + str(len(joint_idx)))
    
    if filter_obs:
        change_tab = cluster_hypergeom_test(aver_all[0].loc[joint_idx], clusters, psi_bin)
        change_tab_reads = cluster_hypergeom_test(aver_all_reads[0].loc[joint_idx], clusters, psi_bin)
    else:
        change_tab = cluster_hypergeom_test(PSI_tab.loc[joint_idx], clusters, psi_bin)
        
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

def test_exon_bimodal_spearman(PSI_tab, exon, obs_min, pca_clust, clusters = 'AC', one_test = False):
    
#     print(psi_lim)
    
    obs_cells = PSI_tab.loc[exon].dropna().index
    
    if one_test:
        
        
        psi = list(PSI_tab.loc[exon, obs_cells])
        pseudotime =  list(pca_clust.loc[obs_cells, 'pseudotime'])
        
        spearman_p = spearmanr(pseudotime, psi)[1]
    
            
    else:
        
        pvals = []
#         print('len')
#         print(len(pca_clust[clusters].unique()))
        for i in range(len(pca_clust[clusters].unique()) - 1):
            
            c1 = i
            c2 = i+1
#             t = time.time()
            c1_ = pca_clust.loc[pca_clust[clusters] == c1].index
            c2_ = pca_clust.loc[pca_clust[clusters] == c2].index
            c1_cells = [x for x in obs_cells if x in c1_]
            c2_cells = [x for x in obs_cells if x in c2_]
#             print('time')
#             print(time.time()-t)
            
            cluster_cells = list(c1_cells) + list(c2_cells)
            psi = list(PSI_tab.loc[exon, cluster_cells])
            pseudotime =  list(pca_clust.loc[cluster_cells, 'pseudotime'])
            
            if ((len(c1_cells)/len(c1_)) >= obs_min) and ((len(c1_cells)/len(c1_)) >= obs_min):
                pvals.append(spearmanr(pseudotime, psi)[1])
                
                
                
#         print(len(pvals))
        if len(pvals) >= 2:
            spearman_p = combine_pvalues(pvals)[1]
        else:
            spearman_p = np.nan
        
    return spearman_p, np.mean(pseudotime), np.mean(psi)





def cluster_spearman_test(PSI_tab, pca_clust, clusters, obs_min, correction = 'fdr_bh', 
                          correct_multitest = True, one_test = False):
    cluster_array = []
    pvals = []
    exon_pass = []
    not_pass = 0
    for i in tqdm(range(len(PSI_tab.index))):

        exon = PSI_tab.index[i]
        chi_p, pos, neg = test_exon_bimodal_spearman(PSI_tab, exon, obs_min, pca_clust, clusters = 'AC', 
                                                     one_test = one_test)
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
#     if correct_multitest:
#         pvals_adj = multipletests(pvals, method=correction)[1]
#     else:
#         pvals_adj = pvals
    cluster_df = pd.DataFrame()
    cluster_df['pvals'] = pvals
    cluster_df.index = exon_pass
    
    return cluster_df


def test_spearman_filters(PSI_tab, mrna_counts, mrna_per_event, read_counts, coverage_tab, pca_clust, clusters = 'AC',
                     psi_min = 0.2, obs_min = 0.5, mrna_min = 10, read_min = 10, filter_obs = False, 
                    dset_name = '', one_test = False, correct_multitest=False):
    
    observed = PSI_tab.loc[PSI_tab.isna().mean(axis=1) <= (1-obs_min)].index
    int_genes, int_exons = spu.get_int_events(PSI_tab.loc[observed], mrna_counts, psi_min)


    aver = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, mrna_min, 0, obs_min)
    
    aver_mrna_only = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, mrna_min, 0, obs_min, filter_cj=False)
    
    aver_reads = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, 0, read_min, obs_min)
    
    aver_all = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, mrna_min, 0, 0)
    
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
                                                  one_test = one_test)
        change_tab_reads, chi_p_reads = cluster_spearman_test(aver_all_reads[0].loc[joint_idx], pca_clust, 
                                                              clusters, one_test = one_test)
        change_tab_reads['pvals'] = chi_p_reads
    else:
        change_tab = cluster_spearman_test(PSI_tab.loc[joint_idx], pca_clust, clusters, obs_min, one_test = one_test)
        
    
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
    For a given p-value cutoff x: 
    — The set of “True instances” X1 is all exons with p-value <= x
    — The set of “False instances” X2 is all exons with p-value > x
    — The set of “Positive calls”  X3 is all exons that you keep after filtering
    — The set of “Negative calls”  X4 is all exons that you throw out after filtering

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

    recall = len(TP)/len(X1) # Also called sensitivity, or TPR
    precision = len(TP)/len(X3)
    specificity = len(TN)/len(X2)
    F1_score = (1+(beta**2))*(precision*recall)/(((beta**2)*precision)+recall)
    hyper_recall_specificity = (1+(beta**2))*(specificity*recall)/(((beta**2)*specificity)+recall)
    
    
    LRplus = recall/(len(FP)/len(X2))
    LRminus = (len(FN)/len(X1))/(specificity)
    try:
        DOR = (recall*specificity)/((1-recall)*(1-specificity))
    except:
        DOR = np.nan
        
    ACC = (len(TP) + len(TN)) / (len(X1)+len(X2))
    BA = (recall*specificity)/2
    
    return recall, precision, specificity, F1_score, hyper_recall_specificity, LRplus, LRminus, DOR, ACC, BA


def summary_curves(pvals, selected_exons, p_limits, beta):
    recall_list = []
    precision_list = []
    specificity_list = []
    f1_score_list = []
    hrs_list = []
    LRplus_list = []
    LRminus_list = []
    DOR_list = []
    ACC_list = []
    BA_list = []
    
    for x in p_limits:
#     for i in tqdm(range(len(p_limits))):
#         x = p_limits[i]
        
        recall, precision, specificity, f1_score, hrs, LRplus, LRminus, DOR , ACC, BA = get_statistics(pvals, selected_exons, x, beta)
        recall_list.append(recall)
        precision_list.append(precision)
        specificity_list.append(specificity)
        f1_score_list.append(f1_score)
        hrs_list.append(hrs)
        LRplus_list.append(LRplus)
        LRminus_list.append(LRminus)
#         DOR = LRplus/LRminus
        DOR_list.append(DOR)
        ACC_list.append(ACC)
        BA_list.append(BA)
        
    return recall_list, precision_list, specificity_list, f1_score_list, hrs_list, LRplus_list, LRminus_list, DOR_list, ACC_list, BA_list
    


def summary_plots(pvals, selected_mrna, selected_mrna_only, selected_read, dset_name, 
                  p_lim=1, p_low = 0.0001, p_steps = 0.0001, linear_p = True, beta=1, plot_line=True):

    if linear_p:
        p_limits = np.arange(p_low, p_lim + p_steps, p_steps)
        
    else:
        p_limits = 10**(-np.arange(-np.log10(p_lim), -np.log10(p_low)+p_steps, p_steps))
               #10**np.arange(-np.log10(p_lim), -np.log10(p_low)+p_steps, p_steps)
    
    recall_mrna, precision_mrna, specificity_mrna, f1_score_mrna, hrs_mrna, LRplus_mrna, LRminus_mrna, DOR_mrna, ACC_mrna, BA_mrna = summary_curves(pvals, selected_mrna, p_limits, beta)
#     recall_mrna_only, precision_mrna_only, specificity_mrna_only, f1_score_mrna_only = summary_curves(pvals, selected_mrna_only, 
#                                                                                                       p_limits, beta)
    recall_read, precision_read, specificity_read, f1_score_read, hrs_read, LRplus_read, LRminus_read, DOR_read, ACC_read, BA_read = summary_curves(pvals, selected_read, p_limits, beta)
    
    '''
    plt.figure()
    plt.scatter(-np.log10(p_limits), recall_read, c='darkred', label='read filter')
    plt.scatter(-np.log10(p_limits), recall_mrna_only, c='forestgreen', label='mrna only filter')
    plt.scatter(-np.log10(p_limits), recall_mrna, c='navy', label='mrna filter')
    
    plt.title('Recall curves, ' + dset_name)
    plt.xlabel('-log10 p-value')
    plt.ylabel('Recall')
    plt.legend(frameon=False)
    plt.show()
    
    plt.figure()
    plt.scatter(-np.log10(p_limits), precision_read, c='darkred', label='read filter')
    plt.scatter(-np.log10(p_limits), precision_mrna_only, c='forestgreen', label='mrna only filter')
    plt.scatter(-np.log10(p_limits), precision_mrna, c='navy', label='mrna filter')
    
    plt.title('Precision curves, ' + dset_name)
    plt.xlabel('-log10 p-value')
    plt.ylabel('Precision')
    plt.legend(frameon=False)
    plt.show()
    
    plt.figure()
    plt.scatter(-np.log10(p_limits), specificity_read, c='darkred', label='read filter')
    plt.scatter(-np.log10(p_limits), specificity_mrna_only, c='forestgreen', label='mrna only filter')
    plt.scatter(-np.log10(p_limits), specificity_mrna, c='navy', label='mrna filter')
    
    plt.title('Specificity curves, ' + dset_name)
    plt.xlabel('-log10 p-value')
    plt.ylabel('Specificity')
    plt.legend(frameon=False)
    plt.show()
    
    '''
    figsize(14,10)
    fig = plt.figure()
    ax  = plt.subplot(1,1,1)
    if plot_line:
    
        ax.plot(-np.log10(p_limits), f1_score_read, c='darkred', label='read filter', linewidth=5)
    #     plt.scatter(-np.log10(p_limits), f1_score_mrna_only, c='forestgreen', label='mrna only filter')
        ax.plot(-np.log10(p_limits), f1_score_mrna, c='navy', label='mrna filter', linewidth=5)
    
    else:
        ax.scatter(-np.log10(p_limits), f1_score_read, c='darkred', label='read filter',s=250, edgecolors='none')
    #     plt.scatter(-np.log10(p_limits), f1_score_mrna_only, c='forestgreen', label='mrna only filter')
        ax.scatter(-np.log10(p_limits), f1_score_mrna, c='navy', label='mrna filter', s=250, edgecolors='none')
    
        
    
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(labelsize=42, length=5)
    plt.title('F'+str(round(beta, 2))+' curves, ' + dset_name, fontsize=42)
    plt.xlabel('-log10 p-value', fontsize=42)
    plt.ylabel('F'+str(round(beta, 2))+' score', fontsize=42)
    plt.legend(frameon=False, fontsize=42)
    
    plt.show()
    
    
    figsize(14,10)
    fig = plt.figure()
    ax  = plt.subplot(1,1,1)
    if plot_line:
    
        ax.plot(-np.log10(p_limits), hrs_read, c='darkred', label='read filter', linewidth=5)
    #     plt.scatter(-np.log10(p_limits), f1_score_mrna_only, c='forestgreen', label='mrna only filter')
        ax.plot(-np.log10(p_limits), hrs_mrna, c='navy', label='mrna filter', linewidth=5)
    
    else:
        ax.scatter(-np.log10(p_limits), hrs_read, c='darkred', label='read filter',s=250, edgecolors='none')
    #     plt.scatter(-np.log10(p_limits), f1_score_mrna_only, c='forestgreen', label='mrna only filter')
        ax.scatter(-np.log10(p_limits), hrs_mrna, c='navy', label='mrna filter', s=250, edgecolors='none')
    
        
    
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(labelsize=42, length=5)
    plt.title('Sensitivity/specificity, ' + dset_name, fontsize=42)
    plt.xlabel('-log10 p-value', fontsize=42)
    plt.ylabel('Harmonic mean', fontsize=42)
    plt.legend(frameon=False, fontsize=42)
    
    plt.show()
    
    
    
    ##################
    
    figsize(14,10)
    fig = plt.figure()
    ax  = plt.subplot(1,1,1)
    if plot_line:
    
        ax.plot(-np.log10(p_limits), LRplus_read, c='darkred', label='read filter', linewidth=5)
    #     plt.scatter(-np.log10(p_limits), f1_score_mrna_only, c='forestgreen', label='mrna only filter')
        ax.plot(-np.log10(p_limits), LRplus_mrna, c='navy', label='mrna filter', linewidth=5)
    
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(labelsize=42, length=5)
    plt.title('Positive likelihood ratio, ' + dset_name, fontsize=42)
    plt.xlabel('-log10 p-value', fontsize=42)
    plt.ylabel('LR+', fontsize=42)
    plt.legend(frameon=False, fontsize=42)
    
    plt.show()

    
    
    figsize(14,10)
    fig = plt.figure()
    ax  = plt.subplot(1,1,1)
    if plot_line:
    
        ax.plot(-np.log10(p_limits), LRminus_read, c='darkred', label='read filter', linewidth=5)
    #     plt.scatter(-np.log10(p_limits), f1_score_mrna_only, c='forestgreen', label='mrna only filter')
        ax.plot(-np.log10(p_limits), LRminus_mrna, c='navy', label='mrna filter', linewidth=5)
    
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(labelsize=42, length=5)
    plt.title('Negative likelihood ratio, ' + dset_name, fontsize=42)
    plt.xlabel('-log10 p-value', fontsize=42)
    plt.ylabel('LR-', fontsize=42)
    plt.legend(frameon=False, fontsize=42)
    
    plt.show()
    
    
    
    figsize(14,10)
    fig = plt.figure()
    ax  = plt.subplot(1,1,1)
    if plot_line:
    
        ax.plot(-np.log10(p_limits), DOR_read, c='darkred', label='read filter', linewidth=5)
    #     plt.scatter(-np.log10(p_limits), f1_score_mrna_only, c='forestgreen', label='mrna only filter')
        ax.plot(-np.log10(p_limits), DOR_mrna, c='navy', label='mrna filter', linewidth=5)
    
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(labelsize=42, length=5)
    plt.title('Diagnostic odds ratio, ' + dset_name, fontsize=42)
    plt.xlabel('-log10 p-value', fontsize=42)
    plt.ylabel('DOR', fontsize=42)
    plt.legend(frameon=False, fontsize=42)
    
    plt.show()
    
    
    
    
    
    
    figsize(14,10)
    fig = plt.figure()
    ax  = plt.subplot(1,1,1)
    if plot_line:
    
        ax.plot(-np.log10(p_limits), ACC_read, c='darkred', label='read filter', linewidth=5)
    #     plt.scatter(-np.log10(p_limits), f1_score_mrna_only, c='forestgreen', label='mrna only filter')
        ax.plot(-np.log10(p_limits), ACC_mrna, c='navy', label='mrna filter', linewidth=5)
    
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(labelsize=42, length=5)
    plt.title('Accuracy, ' + dset_name, fontsize=42)
    plt.xlabel('-log10 p-value', fontsize=42)
    plt.ylabel('ACC', fontsize=42)
    plt.legend(frameon=False, fontsize=42)
    
    plt.show()
    
    
    
    
    
    figsize(14,10)
    fig = plt.figure()
    ax  = plt.subplot(1,1,1)
    if plot_line:
    
        ax.plot(-np.log10(p_limits), BA_read, c='darkred', label='read filter', linewidth=5)
    #     plt.scatter(-np.log10(p_limits), f1_score_mrna_only, c='forestgreen', label='mrna only filter')
        ax.plot(-np.log10(p_limits), BA_mrna, c='navy', label='mrna filter', linewidth=5)
    
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(labelsize=42, length=5)
    plt.title('Balanced Accuracy, ' + dset_name, fontsize=42)
    plt.xlabel('-log10 p-value', fontsize=42)
    plt.ylabel('BA', fontsize=42)
    plt.legend(frameon=False, fontsize=42)
    
    plt.show()


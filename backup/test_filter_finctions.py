from scipy.stats import probplot

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

        clust_cells = pca_clust.loc[pca_clust[clusters] == i].index
        c_cells = [x for x in obs_cells if x in clust_cells]
        
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
    
    for i in tqdm(range(len(PSI_tab.index))):
        
        exon = PSI_tab.index[i]
        anv_p, pos, neg = test_exon_bimodal_anova(PSI_tab, exon, pca_clust, clusters = 'AC', obs_min=obs_min)
        if not np.isnan(anv_p):
            cluster_array.append(pos/(neg+pos))
            pvals.append(anv_p)
            exon_pass.append(exon)

        else:
            not_pass += 1
        
    print('not pass: '+str(not_pass))
    print('tested exons: ' + str(len(pvals)))

    cluster_df = pd.DataFrame()
    cluster_df['pvals'] = pvals
    cluster_df.index = exon_pass
    
    return cluster_df

import random

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
    
    
    
    
    change_tab = cluster_anova_test(aver[0], pca_clust, clusters, obs_min=obs_min)
    change_tab_mrna_only = cluster_anova_test(aver_mrna_only[0], pca_clust, clusters, obs_min=obs_min)
    change_tab_reads = cluster_anova_test(aver_reads[0], pca_clust, clusters, obs_min=obs_min)
    change_tab_unfiltered = cluster_anova_test(PSI_tab.loc[int_exons], pca_clust, clusters, obs_min=obs_min)

#     joint_idx = [x for x in int_exons if x in change_tab.index]
#     joint_idx = [x for x in joint_idx if x in change_tab_mrna_only.index]
#     joint_idx = [x for x in joint_idx if x in change_tab_reads.index]
#     joint_idx = [x for x in joint_idx if x in change_tab_unfiltered.index]
    
    
    missing = PSI_tab.loc[aver[0].index].isna().sum(axis=1)
    missing_mrna = aver[0].isna().sum(axis=1)
#     missing_mrna_only = aver_mrna_only[0].loc[joint_idx].isna().sum(axis=1)
#     missing_reads = aver_reads[0].loc[joint_idx].isna().sum(axis=1)
    
    missing_dif = pd.DataFrame()
    missing_dif['mrna'] = missing_mrna - missing
#     missing_dif['mrna_only'] = missing_mrna_only - missing
#     missing_dif['mrna_read'] = missing_reads - missing
    
    
#     print(missing_dif.max(axis=1))
    
#     aver_random = PSI_tab.loc[joint_idx].copy()
    
    mask = pd.DataFrame()
    
    for exon in aver[0].index:
        non_na_cells = list(PSI_tab.loc[exon].dropna().index)
        remove_n = missing_dif.loc[exon].max()
        remove_cells = r.choice(non_na_cells, int(remove_n))
        filter_col = [x in remove_cells for x in aver[0].columns]
        mask[exon] = filter_col
        
    mask.index = aver[0].columns
    
    aver_random = PSI_tab.loc[aver[0].index].copy()
    aver_random = aver_random.mask(mask.T)
    
    
    change_tab_random = cluster_anova_test(aver_random, pca_clust, clusters, obs_min=obs_min)
    
    
    joint_idx = [x for x in aver[0].index if x in change_tab_random.index]
    
    
    
    
    
    
    
    
    
    print('Selected intermediate exons: ' + str(len(joint_idx)))
    
    figsize(6, 4)
    
#     plt.scatter(probplot(change_tab.loc[joint_idx, 'pvals'], dist='uniform')[0][0], 
#                      probplot(change_tab.loc[joint_idx, 'pvals'], dist='uniform')[0][1], 
#                      alpha = 0.5, c='navy', label='filter')
    
#     plt.scatter(probplot(change_tab_reads.loc[joint_idx, 'pvals'], dist='uniform')[0][0], 
#                      probplot(change_tab_reads.loc[joint_idx, 'pvals'], dist='uniform')[0][1], 
#                      alpha = 0.5, c='darkred', label='read filter')
    
#     plt.scatter(probplot(change_tab_mrna_only.loc[joint_idx, 'pvals'], dist='uniform')[0][0], 
#                      probplot(change_tab_mrna_only.loc[joint_idx, 'pvals'], dist='uniform')[0][1], 
#                      alpha = 0.5, c='purple', label='mrna only')

#     plt.scatter(probplot(change_tab_unfiltered.loc[joint_idx, 'pvals'], dist='uniform')[0][0], 
#                      probplot(change_tab_unfiltered.loc[joint_idx, 'pvals'], dist='uniform')[0][1], 
#                      alpha = 0.5, c='purple', label='mrna only')
    
#     plt.scatter(probplot(change_tab_random.loc[joint_idx, 'pvals'], dist='uniform')[0][0], 
#                      probplot(change_tab_random.loc[joint_idx, 'pvals'], dist='uniform')[0][1], 
#                      alpha = 0.5, c='skyblue', label='mrna only')
    
#     plt.show()
    
#     print((change_tab.loc[joint_idx, 'pvals']<= 0.05).sum()/len(joint_idx))
#     print((change_tab_reads.loc[joint_idx, 'pvals']<= 0.05).sum()/len(joint_idx))
#     print((change_tab_mrna_only.loc[joint_idx, 'pvals']<= 0.05).sum()/len(joint_idx))
#     print((change_tab_unfiltered.loc[joint_idx, 'pvals']<= 0.05).sum()/len(joint_idx))
#     print((change_tab_random.loc[joint_idx, 'pvals']<= 0.05).sum()/len(joint_idx))

    plt.scatter(probplot(change_tab.pvals, dist='uniform')[0][0], 
                     probplot(change_tab.pvals, dist='uniform')[0][1], 
                     alpha = 0.5, c='navy', label='filter')
    
    plt.scatter(probplot(change_tab_reads.pvals, dist='uniform')[0][0], 
                     probplot(change_tab_reads.pvals, dist='uniform')[0][1], 
                     alpha = 0.5, c='darkred', label='read filter')
    
    plt.scatter(probplot(change_tab_mrna_only.pvals, dist='uniform')[0][0], 
                     probplot(change_tab_mrna_only.pvals, dist='uniform')[0][1], 
                     alpha = 0.5, c='purple', label='mrna only')
    
    plt.scatter(probplot(change_tab_unfiltered.pvals, dist='uniform')[0][0], 
                     probplot(change_tab_unfiltered.pvals, dist='uniform')[0][1], 
                     alpha = 0.5, c='darkorange', label='unfiltered')
    
    
    plt.scatter(probplot(change_tab_random.pvals, dist='uniform')[0][0], 
                     probplot(change_tab_random.pvals, dist='uniform')[0][1], 
                     alpha = 0.5, c='skyblue', label='unfiltered')
    
    plt.show()
    
    print((change_tab.pvals<= 0.05).mean())
    print((change_tab_reads.pvals<= 0.05).mean())
    print((change_tab_mrna_only.pvals<= 0.05).mean())
    print((change_tab_unfiltered.pvals<= 0.05).mean())
    print((change_tab_random.pvals<= 0.05).mean())

    plt.hist(change_tab.pvals - change_tab_unfiltered.loc[change_tab.index, 'pvals'], 
             color='navy', alpha=0.5, density=True, bins=20)
#     plt.hist(change_tab_reads.pvals - change_tab_unfiltered.loc[change_tab_reads.index, 'pvals'], 
#              color='darkred', alpha=0.5, density=True, bins=20)
#     plt.hist(change_tab_mrna_only.pvals - change_tab_unfiltered.loc[change_tab_mrna_only.index, 'pvals'], 
#              color='forestgreen', alpha=0.5, density=True, bins=20)
#     plt.hist(change_tab_random.pvals - change_tab_unfiltered.loc[change_tab_random.index, 'pvals'], 
#              color='purple', alpha=0.5, density=True, bins=20)

    plt.show()
############################





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
    
    
    change_tab = cluster_chi_test(aver[0], clusters, psi_bin)
    change_tab_mrna_only = cluster_chi_test(aver_mrna_only[0], clusters, psi_bin)
    change_tab_reads = cluster_chi_test(aver_reads[0], clusters, psi_bin)
    change_tab_unfiltered = cluster_chi_test(PSI_tab.loc[int_exons], clusters, psi_bin)

    joint_idx = [x for x in int_exons if x in change_tab.index]
    joint_idx = [x for x in joint_idx if x in change_tab_mrna_only.index]
    joint_idx = [x for x in joint_idx if x in change_tab_reads.index]
    joint_idx = [x for x in joint_idx if x in change_tab_unfiltered.index]
    
    
    
    joint_idx = [x for x in int_exons if x in change_tab.index]
    joint_idx = [x for x in joint_idx if x in change_tab_mrna_only.index]
    joint_idx = [x for x in joint_idx if x in change_tab_reads.index]
    joint_idx = [x for x in joint_idx if x in change_tab_unfiltered.index]
    
    
    missing = PSI_tab.loc[joint_idx].isna().sum(axis=1)
    missing_mrna = aver[0].loc[joint_idx].isna().sum(axis=1)
    missing_mrna_only = aver_mrna_only[0].loc[joint_idx].isna().sum(axis=1)
    missing_reads = aver_reads[0].loc[joint_idx].isna().sum(axis=1)
    
    missing_dif = pd.DataFrame()
    missing_dif['mrna'] = missing_mrna - missing
    missing_dif['mrna_only'] = missing_mrna_only - missing
    missing_dif['mrna_read'] = missing_reads - missing
    
    
#     print(missing_dif.max(axis=1))
    
#     aver_random = PSI_tab.loc[joint_idx].copy()
    
    mask = pd.DataFrame()
    
    for exon in joint_idx:
        non_na_cells = list(PSI_tab.loc[exon].dropna().index)
        remove_n = missing_dif.loc[exon].max()
        remove_cells = r.choice(non_na_cells, int(remove_n))
        filter_col = [x in remove_cells for x in aver[0].columns]
        mask[exon] = filter_col
        
    mask.index = aver[0].columns
    
    aver_random = PSI_tab.loc[joint_idx].copy()
    aver_random = aver_random.mask(mask.T)
    
    
    change_tab_random = cluster_chi_test(aver_random, clusters, psi_bin)
    
    
    joint_idx = [x for x in joint_idx if x in change_tab_random.index]
    
    
    
    
    
    
    
    
    
    print('Selected intermediate exons: ' + str(len(joint_idx)))
    
    figsize(6, 4)
    
    plt.scatter(probplot(change_tab.loc[joint_idx, 'pvals'], dist='uniform')[0][0], 
                     probplot(change_tab.loc[joint_idx, 'pvals'], dist='uniform')[0][1], 
                     alpha = 0.5, c='navy', label='filter')
    
    plt.scatter(probplot(change_tab_reads.loc[joint_idx, 'pvals'], dist='uniform')[0][0], 
                     probplot(change_tab_reads.loc[joint_idx, 'pvals'], dist='uniform')[0][1], 
                     alpha = 0.5, c='darkred', label='read filter')
    
    plt.scatter(probplot(change_tab_mrna_only.loc[joint_idx, 'pvals'], dist='uniform')[0][0], 
                     probplot(change_tab_mrna_only.loc[joint_idx, 'pvals'], dist='uniform')[0][1], 
                     alpha = 0.5, c='purple', label='mrna only')

    plt.scatter(probplot(change_tab_unfiltered.loc[joint_idx, 'pvals'], dist='uniform')[0][0], 
                     probplot(change_tab_unfiltered.loc[joint_idx, 'pvals'], dist='uniform')[0][1], 
                     alpha = 0.5, c='purple', label='mrna only')
    
    plt.scatter(probplot(change_tab_random.loc[joint_idx, 'pvals'], dist='uniform')[0][0], 
                     probplot(change_tab_random.loc[joint_idx, 'pvals'], dist='uniform')[0][1], 
                     alpha = 0.5, c='skyblue', label='mrna only')
    
    plt.show()
    
    print((change_tab.loc[joint_idx, 'pvals']<= 0.05).sum()/len(joint_idx))
    print((change_tab_reads.loc[joint_idx, 'pvals']<= 0.05).sum()/len(joint_idx))
    print((change_tab_mrna_only.loc[joint_idx, 'pvals']<= 0.05).sum()/len(joint_idx))
    print((change_tab_unfiltered.loc[joint_idx, 'pvals']<= 0.05).sum()/len(joint_idx))
    print((change_tab_random.loc[joint_idx, 'pvals']<= 0.05).sum()/len(joint_idx))

#     plt.scatter(probplot(change_tab.pvals, dist='uniform')[0][0], 
#                      probplot(change_tab.pvals, dist='uniform')[0][1], 
#                      alpha = 0.5, c='navy', label='filter')
    
#     plt.scatter(probplot(change_tab_reads.pvals, dist='uniform')[0][0], 
#                      probplot(change_tab_reads.pvals, dist='uniform')[0][1], 
#                      alpha = 0.5, c='darkred', label='read filter')
    
#     plt.scatter(probplot(change_tab_mrna_only.pvals, dist='uniform')[0][0], 
#                      probplot(change_tab_mrna_only.pvals, dist='uniform')[0][1], 
#                      alpha = 0.5, c='purple', label='mrna only')
    
#     plt.scatter(probplot(change_tab_unfiltered.pvals, dist='uniform')[0][0], 
#                      probplot(change_tab_unfiltered.pvals, dist='uniform')[0][1], 
#                      alpha = 0.5, c='darkorange', label='unfiltered')
    
#     plt.show()
    
#     print((change_tab.pvals<= 0.05).mean())
#     print((change_tab_reads.pvals<= 0.05).mean())
#     print((change_tab_mrna_only.pvals<= 0.05).mean())
#     print((change_tab_unfiltered.pvals<= 0.05).mean())


def test_spearman_filters(PSI_tab, mrna_counts, mrna_per_event, read_counts, coverage_tab, pca_clust, clusters = 'AC',
                     psi_min = 0.2, obs_min = 0.5, mrna_min = 10, read_min = 10, filter_obs = False, 
                    dset_name = '', one_test = False, correct_multitest=False):
    
    
    
    observed = PSI_tab.loc[PSI_tab.isna().mean(axis=1) <= (1-obs_min)].index
    int_genes, int_exons = spu.get_int_events(PSI_tab.loc[observed], mrna_counts, psi_min)


    aver = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, mrna_min, 0, obs_min)
    
    aver_mrna_only = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, mrna_min, 0, obs_min, filter_cj = False)
    
    aver_reads = filter_psi(PSI_tab, int_exons, mrna_per_event, coverage_tab['SJ_coverage'],
              read_counts, 0, read_min, obs_min)
    
    
    change_tab = cluster_spearman_test(aver[0], pca_clust, clusters, obs_min, one_test = one_test)
    change_tab_mrna_only = cluster_spearman_test(aver_mrna_only[0], pca_clust, clusters, obs_min, one_test = one_test)
    change_tab_reads = cluster_spearman_test(aver_reads[0], pca_clust, clusters, obs_min, one_test = one_test)
    change_tab_unfiltered = cluster_spearman_test(PSI_tab.loc[int_exons], pca_clust, clusters, obs_min, one_test = one_test)

    joint_idx = [x for x in int_exons if x in change_tab.index]
    joint_idx = [x for x in joint_idx if x in change_tab_mrna_only.index]
    joint_idx = [x for x in joint_idx if x in change_tab_reads.index]
    joint_idx = [x for x in joint_idx if x in change_tab_unfiltered.index]
    
    
    
    joint_idx = [x for x in int_exons if x in change_tab.index]
    joint_idx = [x for x in joint_idx if x in change_tab_mrna_only.index]
    joint_idx = [x for x in joint_idx if x in change_tab_reads.index]
    joint_idx = [x for x in joint_idx if x in change_tab_unfiltered.index]
    
    
    missing = PSI_tab.loc[joint_idx].isna().sum(axis=1)
    missing_mrna = aver[0].loc[joint_idx].isna().sum(axis=1)
    missing_mrna_only = aver_mrna_only[0].loc[joint_idx].isna().sum(axis=1)
    missing_reads = aver_reads[0].loc[joint_idx].isna().sum(axis=1)
    
    missing_dif = pd.DataFrame()
    missing_dif['mrna'] = missing_mrna - missing
    missing_dif['mrna_only'] = missing_mrna_only - missing
    missing_dif['mrna_read'] = missing_reads - missing
    
    
#     print(missing_dif.max(axis=1))
    
#     aver_random = PSI_tab.loc[joint_idx].copy()
    
    mask = pd.DataFrame()
    
    for exon in joint_idx:
        non_na_cells = list(PSI_tab.loc[exon].dropna().index)
        remove_n = missing_dif.loc[exon].max()
        remove_cells = r.choice(non_na_cells, int(remove_n))
        filter_col = [x in remove_cells for x in aver[0].columns]
        mask[exon] = filter_col
        
    mask.index = aver[0].columns
    
    aver_random = PSI_tab.loc[joint_idx].copy()
    aver_random = aver_random.mask(mask.T)
    
    
    change_tab_random = cluster_spearman_test(aver_random, pca_clust, clusters, obs_min, one_test = one_test)
    
    joint_idx = [x for x in joint_idx if x in change_tab_random.index]
    
    
    
    
    
    
    
    
    
    print('Selected intermediate exons: ' + str(len(joint_idx)))
    
    figsize(6, 4)
    
    plt.scatter(probplot(change_tab.loc[joint_idx, 'pvals'], dist='uniform')[0][0], 
                     probplot(change_tab.loc[joint_idx, 'pvals'], dist='uniform')[0][1], 
                     alpha = 0.5, c='navy', label='filter')
    
    plt.scatter(probplot(change_tab_reads.loc[joint_idx, 'pvals'], dist='uniform')[0][0], 
                     probplot(change_tab_reads.loc[joint_idx, 'pvals'], dist='uniform')[0][1], 
                     alpha = 0.5, c='darkred', label='read filter')
    
    plt.scatter(probplot(change_tab_mrna_only.loc[joint_idx, 'pvals'], dist='uniform')[0][0], 
                     probplot(change_tab_mrna_only.loc[joint_idx, 'pvals'], dist='uniform')[0][1], 
                     alpha = 0.5, c='forestgreen', label='mrna only')

    plt.scatter(probplot(change_tab_unfiltered.loc[joint_idx, 'pvals'], dist='uniform')[0][0], 
                     probplot(change_tab_unfiltered.loc[joint_idx, 'pvals'], dist='uniform')[0][1], 
                     alpha = 0.5, c='purple', label='mrna only')
    
    plt.scatter(probplot(change_tab_random.loc[joint_idx, 'pvals'], dist='uniform')[0][0], 
                     probplot(change_tab_random.loc[joint_idx, 'pvals'], dist='uniform')[0][1], 
                     alpha = 0.5, c='skyblue', label='mrna only')
    
    plt.show()
    
    print((change_tab.loc[joint_idx, 'pvals']<= 0.05).sum()/len(joint_idx))
    print((change_tab_reads.loc[joint_idx, 'pvals']<= 0.05).sum()/len(joint_idx))
    print((change_tab_mrna_only.loc[joint_idx, 'pvals']<= 0.05).sum()/len(joint_idx))
    print((change_tab_unfiltered.loc[joint_idx, 'pvals']<= 0.05).sum()/len(joint_idx))
    print((change_tab_random.loc[joint_idx, 'pvals']<= 0.05).sum()/len(joint_idx))
    
    
    plt.hist(change_tab.pvals - change_tab_unfiltered.loc[change_tab.index, 'pvals'], 
             color='navy', alpha=0.5, density=True, bins=20)
    plt.hist(change_tab_reads.pvals - change_tab_unfiltered.loc[change_tab_reads.index, 'pvals'], 
             color='darkred', alpha=0.5, density=True, bins=20)
    plt.hist(change_tab_mrna_only.pvals - change_tab_unfiltered.loc[change_tab_mrna_only.index, 'pvals'], 
             color='forestgreen', alpha=0.5, density=True, bins=20)
    plt.hist(change_tab_random.pvals - change_tab_unfiltered.loc[change_tab_random.index, 'pvals'], 
             color='purple', alpha=0.5, density=True, bins=20)

    plt.show()
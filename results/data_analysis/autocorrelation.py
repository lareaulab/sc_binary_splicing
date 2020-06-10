data_dir = '/mnt/lareaulab/cfbuenabadn/SingleCell/data/'

######################
# load_data_short.py #
######################

print('loading data')

import numpy as np
import pandas as pd
import os
import matplotlib.cm as cm
from matplotlib import pyplot as plt
from scipy import stats as st
import seaborn as sns

import numpy.random as r
import sys
sys.path.insert(0, '/mnt/lareaulab/cfbuenabadn/sc_binary_splicing/utils/')
import splicing_utils as spu
from splicing_utils import *
import single_cell_plots as scp
from single_cell_plots import *

import numpy as np
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering

data_dir = '/mnt/lareaulab/cfbuenabadn/SingleCell/data/'

# load PSI tables
chen_PSI = pd.read_csv(data_dir + 'chen/processed_tables/chen.skipped_exons_psi.tab', sep='\t', index_col=0)
song_PSI = pd.read_csv(data_dir + 'song/processed_tables/song.skipped_exons_psi.tab', sep='\t',  index_col=0)
trapnell_PSI = pd.read_csv(data_dir + 'trapnell/processed_tables/trapnell.skipped_exons_psi.tab', 
                    sep='\t',  index_col=0)
lescroart_PSI = pd.read_csv(data_dir + 'lescroart/processed_tables/lescroart.skipped_exons_psi.tab', 
                    sep='\t',  index_col=0)

# SJ read tables
chen_read_counts = pd.read_csv(data_dir + 'chen/processed_tables/chen.skipped_exons_SJreads.tab', sep='\t', index_col=0)
song_read_counts = pd.read_csv(data_dir + 'song/processed_tables/song.skipped_exons_SJreads.tab', sep='\t', index_col=0)
trapnell_read_counts = pd.read_csv(data_dir + 'trapnell/processed_tables/trapnell.skipped_exons_SJreads.tab', 
                            sep='\t', index_col=0)
lescroart_read_counts = pd.read_csv(data_dir + 'lescroart/processed_tables/lescroart.skipped_exons_SJreads.tab', 
                             sep='\t', index_col=0)

# mRNA tables
chen_mrna_counts = pd.read_csv(data_dir + 'chen/processed_tables/chen.mrna_counts.tab', sep='\t', index_col=0)
song_mrna_counts = pd.read_csv(data_dir + 'song/processed_tables/song.mrna_counts.tab', sep='\t', index_col=0)
trapnell_mrna_counts = pd.read_csv(data_dir + 'trapnell/processed_tables/trapnell.mrna_counts.tab', 
                            sep='\t', index_col=0)
lescroart_mrna_counts = pd.read_csv(data_dir + 'lescroart/processed_tables/lescroart.mrna_counts.tab', 
                             sep='\t', index_col=0)

# mRNA per event rables
mrna_per_event_chen = pd.read_csv(data_dir + 'chen/processed_tables/chen.mrna_counts_per_event.tab', sep='\t', index_col=0)
mrna_per_event_song = pd.read_csv(data_dir + 'song/processed_tables/song.mrna_counts_per_event.tab', sep='\t', index_col=0)
mrna_per_event_trapnell = pd.read_csv(data_dir + 'trapnell/processed_tables/trapnell.mrna_counts_per_event.tab', 
                            sep='\t', index_col=0)
mrna_per_event_lescroart = pd.read_csv(data_dir + 'lescroart/processed_tables/lescroart.mrna_counts_per_event.tab', 
                             sep='\t', index_col=0)

# read coverage tables
chen_coverage_tab = pd.read_csv(data_dir + 'chen/processed_tables/chen.read_coverage.tab', 
                          sep='\t', index_col=0)
song_coverage_tab = pd.read_csv(data_dir + 'song/processed_tables/song.read_coverage.tab', 
                          sep='\t', index_col=0)
trapnell_coverage_tab = pd.read_csv(data_dir + 'trapnell/processed_tables/trapnell.read_coverage.tab', 
                          sep='\t', index_col=0)
lescroart_coverage_tab = pd.read_csv(data_dir + 'lescroart/processed_tables/lescroart.read_coverage.tab', 
                          sep='\t', index_col=0)


chen_pca = pd.read_csv(data_dir + 'chen/chen.pca.tab', sep='\t', index_col=0)
chen_pca = chen_pca.sort_values('pseudotime')
chen_pca.PC2 = chen_pca.PC2
chen_pca.line_2 = chen_pca.line_2
chen_index = [x for x in chen_pca.sort_values('pseudotime').index if x in mrna_per_event_chen.columns]

chen_pca = chen_pca.loc[chen_index]
chen_PSI = chen_PSI[chen_index]
mrna_per_event_chen = mrna_per_event_chen[chen_index]
chen_read_counts = chen_read_counts[chen_index]
chen_coverage_tab = chen_coverage_tab.loc[chen_index]

chen_ES2i = chen_pca.loc[chen_pca.cell_type == 'ES2i'].index
chen_ES = chen_pca.loc[chen_pca.cell_type == 'ES'].index
chen_Epi = chen_pca.loc[chen_pca.cell_type == 'Epi'].index
chen_MN = chen_pca.loc[chen_pca.cell_type == 'Motor neuron'].index

###

song_pca = pd.read_csv(data_dir + 'song/song.pca_top_fano.tab', sep='\t', index_col=0)
song_pca = song_pca.sort_values('pseudotime')
song_pca.PC2 = -song_pca.PC2
song_index = [x for x in song_pca.sort_values('pseudotime').index if x in mrna_per_event_song.columns]

song_pca = song_pca.loc[song_index]
song_PSI = song_PSI[song_index]
mrna_per_event_song = mrna_per_event_song[song_index]
song_read_counts = song_read_counts[song_index]
song_coverage_tab = song_coverage_tab.loc[song_index]

song_iPSC = song_pca.loc[song_pca.cell_type == 'iPSC'].index
song_NPC = song_pca.loc[song_pca.cell_type == 'NPC'].index
song_MN = song_pca.loc[song_pca.cell_type == 'MN'].index

###

trapnell_pca = pd.read_csv(data_dir + 'trapnell/trapnell.pca.tab', sep='\t', index_col=0)
trapnell_pca = trapnell_pca.sort_values('pseudotime')
trapnell_pca.PC2 = -trapnell_pca.PC2
trapnell_pca.line_2 = -trapnell_pca.line_2
trapnell_index = [x for x in trapnell_pca.sort_values('pseudotime').index if x in mrna_per_event_trapnell.columns]

trapnell_pca = trapnell_pca.loc[trapnell_index]
trapnell_PSI = trapnell_PSI[trapnell_index]
mrna_per_event_trapnell = mrna_per_event_trapnell[trapnell_index]
trapnell_read_counts = trapnell_read_counts[trapnell_index]
trapnell_coverage_tab = trapnell_coverage_tab.loc[trapnell_index]

trapnell_M00 = trapnell_pca.loc[trapnell_pca.cell_type==0].index
trapnell_M24 = trapnell_pca.loc[trapnell_pca.cell_type==24].index
trapnell_M48 = trapnell_pca.loc[trapnell_pca.cell_type==48].index
trapnell_M72 = trapnell_pca.loc[trapnell_pca.cell_type==72].index


###
        
lescroart_pca = pd.read_csv(data_dir + 'lescroart/lescroart.pca_meta.tab', sep='\t', index_col = 0)
lescroart_index = [x for x in lescroart_pca.index if x in mrna_per_event_lescroart.columns]

lescroart_pca = lescroart_pca.loc[lescroart_index]
lescroart_PSI = lescroart_PSI[lescroart_index]
mrna_per_event_lescroart = mrna_per_event_lescroart[lescroart_index]
lescroart_read_counts = lescroart_read_counts[lescroart_index]
lescroart_coverage_tab = lescroart_coverage_tab.loc[lescroart_index]

lescroart_E6 = lescroart_pca.loc[lescroart_pca.cell_type=='E6.75'].index
lescroart_E7 = lescroart_pca.loc[lescroart_pca.cell_type=='E7.25'].index

def process_subpop(subpop, psi, mrna, mrna_per_event, reads, cj, psi_min = 0.2, mrna_min=10, reads_min = 0, cell_min = 0.5, nbins=11,
                  filter_cj = True):

    int_genes, int_exons = spu.get_int_events(psi[subpop], mrna[subpop], psi_min)
    #print(len(int_genes))
    int_exons = [x for x in int_exons if x in mrna_per_event.index]
    PSI_filtered, PSI_mrna_filtered, good_exons, mrna_filtered, reads_filtered = filter_psi(psi[subpop], int_exons, 
                                                                     mrna_per_event[subpop], cj.loc[subpop], 
                                                                     reads[subpop], mrna_min, reads_min = reads_min,
                                                                                            cell_min=cell_min, filter_cj=filter_cj)

    good_cells = PSI_filtered.dropna(axis=1, how='all').columns
    good_subpop = [x for x in subpop if x in good_cells]
    PSI_good = PSI_filtered[good_cells]

    hist_complete, hist_intermediate = scp.get_bins_table2(PSI_filtered[good_subpop], mrna_filtered[good_subpop], nbins)
    hist_complete_exp, hist_intermediate_exp = scp.get_bins_table(PSI_filtered[good_subpop], mrna_filtered[good_subpop])
    
    return PSI_filtered, good_exons, mrna_filtered, reads_filtered, hist_complete, hist_complete_exp
    
# Chen
ac = AgglomerativeClustering(n_clusters=5)
ac_clusters = ac.fit_predict(chen_pca[['PC1', 'PC2']])

chen_pca_clust = chen_pca.copy()
chen_pca_clust['AC'] = ac_clusters

chen_clust_filter = []
for cluster in chen_pca_clust.groupby('AC')['pseudotime'].mean().sort_values().index:
    clust_subpop = chen_pca_clust.index[chen_pca_clust.AC == cluster]
    
    chen_filter = process_subpop(clust_subpop, chen_PSI, chen_mrna_counts, mrna_per_event_chen, 
                                 chen_read_counts, chen_coverage_tab['SJ_coverage'], 0.1, 10, 0, cell_min=0.5)
    
    chen_clust_filter.append(chen_filter)
    
# Song
ac = AgglomerativeClustering(n_clusters=3)
ac_clusters = ac.fit_predict(song_pca[['PC1', 'PC2']])

song_pca_clust = song_pca.copy()
song_pca_clust['AC'] = ac_clusters

song_clust_filter = []
for cluster in song_pca_clust.groupby('AC')['pseudotime'].mean().sort_values().index:
    clust_subpop = song_pca_clust.index[song_pca_clust.AC == cluster]
    
    song_filter = process_subpop(clust_subpop, song_PSI, song_mrna_counts, mrna_per_event_song, 
                                 song_read_counts, song_coverage_tab['SJ_coverage'], 0.1, 10, 0, cell_min=0.5)
    
    song_clust_filter.append(song_filter)
    
    
# Lescroart
ac = AgglomerativeClustering(n_clusters=3)
ac_clusters = ac.fit_predict(lescroart_pca[['PC1', 'PC2']])

lescroart_pca_clust = lescroart_pca.copy()
lescroart_pca_clust['AC'] = ac_clusters


lescroart_clust_filter = []
for cluster in lescroart_pca_clust.cell_type.unique():
    clust_subpop = lescroart_pca_clust.index[lescroart_pca_clust.cell_type == cluster]
    
    lescroart_filter = process_subpop(clust_subpop, lescroart_PSI, lescroart_mrna_counts, mrna_per_event_lescroart, 
                                 lescroart_read_counts, lescroart_coverage_tab['SJ_coverage'], 0.1, 10, 0, cell_min=0.5)
    
    lescroart_clust_filter.append(lescroart_filter)
    
    
###################
# Trapnell
ac = AgglomerativeClustering(n_clusters=4)
ac_clusters = ac.fit_predict(trapnell_pca[['PC1', 'PC2']])

trapnell_pca_clust = trapnell_pca.copy()
trapnell_pca_clust['AC'] = ac_clusters

trapnell_clust_filter = []
for cluster in trapnell_pca_clust.groupby('AC')['pseudotime'].mean().sort_values().index:
    clust_subpop = trapnell_pca_clust.index[trapnell_pca_clust.AC == cluster]
    
    trapnell_filter = process_subpop(clust_subpop, trapnell_PSI, trapnell_mrna_counts, mrna_per_event_trapnell, 
                                 trapnell_read_counts, trapnell_coverage_tab['SJ_coverage'], 0.1, 10, 0, cell_min=0.5)
    
    trapnell_clust_filter.append(trapnell_filter)

#####

from sklearn.decomposition import PCA
from scipy.stats import spearmanr
#import rpy2
#import rpy2.robjects.packages as rpackages
#import rpy2.robjects as robjects
#import rpy2.robjects.numpy2ri as rpyn
from statsmodels.stats.multitest import multipletests
#dt = rpy2.robjects.packages.importr('diptest')

from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering
from scipy.special import logit
from scipy.special import expit
from sklearn.metrics import adjusted_rand_score
from scipy.stats import combine_pvalues


###############

# For the sake of quantification, filter exons between 0.05 and 0.95

###################
# Chen
ac = AgglomerativeClustering(n_clusters=5)
ac_clusters = ac.fit_predict(chen_pca[['PC1', 'PC2']])

# figsize(6,4)
# plt.scatter(chen_pca.PC1, chen_pca.PC2, c=ac_clusters)
# plt.show()

chen_pca_clust = chen_pca.copy()
chen_pca_clust['AC'] = ac_clusters

chen_clust_filter_05 = []
for cluster in chen_pca_clust.groupby('AC')['pseudotime'].mean().sort_values().index:
    clust_subpop = chen_pca_clust.index[chen_pca_clust.AC == cluster]
    
    chen_filter = process_subpop(clust_subpop, chen_PSI, chen_mrna_counts, mrna_per_event_chen, 
                                 chen_read_counts, chen_coverage_tab['SJ_coverage'], 0.05, 10, 0, cell_min=0.5)
    
    chen_clust_filter_05.append(chen_filter)
    
###################
# Song
ac = AgglomerativeClustering(n_clusters=3)
ac_clusters = ac.fit_predict(song_pca[['PC1', 'PC2']])

# figsize(6,4)
# plt.scatter(song_pca.PC1, song_pca.PC2, c=ac_clusters)
# plt.show()

song_pca_clust = song_pca.copy()
song_pca_clust['AC'] = ac_clusters

song_clust_filter_05 = []
for cluster in song_pca_clust.groupby('AC')['pseudotime'].mean().sort_values().index:
    clust_subpop = song_pca_clust.index[song_pca_clust.AC == cluster]
    
    song_filter = process_subpop(clust_subpop, song_PSI, song_mrna_counts, mrna_per_event_song, 
                                 song_read_counts, song_coverage_tab['SJ_coverage'], 0.05, 10, 0, cell_min=0.5)
    
    song_clust_filter_05.append(song_filter)
    
    
###################
# Lescroart
ac = AgglomerativeClustering(n_clusters=3)
ac_clusters = ac.fit_predict(lescroart_pca[['PC1', 'PC2']])

lescroart_clust_filter_05 = []
for cluster in lescroart_pca_clust.cell_type.unique():
    clust_subpop = lescroart_pca_clust.index[lescroart_pca_clust.cell_type == cluster]
    
    lescroart_filter = process_subpop(clust_subpop, lescroart_PSI, lescroart_mrna_counts, mrna_per_event_lescroart, 
                                 lescroart_read_counts, lescroart_coverage_tab['SJ_coverage'], 0.05, 10, 0, cell_min=0.5)
    
    lescroart_clust_filter_05.append(lescroart_filter)
    
    
    
###################
# Trapnell
ac = AgglomerativeClustering(n_clusters=4)
ac_clusters = ac.fit_predict(trapnell_pca[['PC1', 'PC2']])

trapnell_pca_clust = trapnell_pca.copy()
trapnell_pca_clust['AC'] = ac_clusters

trapnell_clust_filter_05 = []
for cluster in trapnell_pca_clust.groupby('AC')['pseudotime'].mean().sort_values().index:
    clust_subpop = trapnell_pca_clust.index[trapnell_pca_clust.AC == cluster]
    
    trapnell_filter = process_subpop(clust_subpop, trapnell_PSI, trapnell_mrna_counts, mrna_per_event_trapnell, 
                                 trapnell_read_counts, trapnell_coverage_tab['SJ_coverage'], 0.05, 10, 0, cell_min=0.5)
    
    trapnell_clust_filter_05.append(trapnell_filter)

    
# %run -i '../../utils/load_data_short.py'

# import sys
# sys.path.insert(0, '../../utils/')
import importlib
# importlib.reload(scp)
# importlib.reload(spu)
# sns.reset_orig()
from scipy.stats import combine_pvalues
# %run -i '../../utils/Kruskal_Wallis_test_functions.py'

from tqdm import tqdm

############

from sklearn.neighbors import NearestNeighbors
from scipy.spatial.distance import pdist
from sklearn.metrics.pairwise import euclidean_distances


print('data loaded')

def get_distance_matrix(pca, k=10):
    nbrs = NearestNeighbors(n_neighbors=k).fit(pca[['PC1', 'PC2']])
    distances, indices = nbrs.kneighbors(pca[['PC1', 'PC2']])
    
    cells = list(pca.index)
    
    W = pd.DataFrame(np.zeros((len(cells), len(cells))))
    W.columns = cells
    W.index = cells
    
    for i in tqdm(range(len(cells))):
        cell_i = cells[i]
        sigma = np.max(distances[i])
        for j in range(len(distances[i])):
            cell_j = cells[indices[i][j]]
            d = distances[i][j]
            w = np.exp(-(d**2)/(sigma**2))        
            W.loc[cell_i, cell_j] = w
    
    return W


def get_signature_matrix(PSI_tab):
    return (PSI_tab - PSI_tab.mean())/PSI_tab.std()


def make_mock_C_scores(norm_PSI, Ws, exon_list, total_cells, mock=100000):
    exon_out_list = []
    C_scores = []
    for i in tqdm(range(mock)):
        mock_run = True
        while mock_run:
            exon = r.choice(exon_list, 1)[0]
#             print(exon)
            scramble_cells = r.choice(norm_PSI.columns, total_cells, replace=False)
            mock_PSI = pd.DataFrame(norm_PSI.loc[exon, scramble_cells]).T

#             print(mock_PSI.shape)
#             print(norm_PSI.shape)

            mock_PSI.columns = norm_PSI.columns
            mock_df = mock_PSI.loc[exon]
#             print(type(mock_df))
#             print(mock_df)
            mock_score = get_C(mock_df, Ws)
            
#             print(mock_score)

            if mock_score >= 0:
                C_scores.append(mock_score)
                exon_out_list.append('mock_'+exon+'_'+str(i))
                mock_run = False
    return exon_out_list, C_scores
                    
def get_C(exon_score, W):
    exon_score = exon_score.dropna()
    obs_cells = exon_score.index
    x = (exon_score.values.reshape(-1, 1) - exon_score.values.reshape(1, -1))
    w = W.loc[obs_cells, obs_cells]
    num = (len(obs_cells)-1)*((w*(x**2)).sum().sum())
    den = (2*w.sum().sum())*np.sum((exon_score - exon_score.mean())**2)
    C = num/den
    score = 1 - C
    return score
    
    
##############################


# Ahora si el bueno

def get_mock_dict(PSI_tab, norm_PSI, Ws, mock=200):

    total_cells = len(PSI_tab.columns)

    exons_05_10 = PSI_tab.index[(np.abs(0.5 - PSI_tab.mean(axis = 1))>0.4) & (np.abs(0.5 - PSI_tab.mean(axis = 1))<=0.45)] 
    exons_10_20 = PSI_tab.index[(np.abs(0.5 - PSI_tab.mean(axis = 1))>0.3) & (np.abs(0.5 - PSI_tab.mean(axis = 1))<=0.40)] 
    exons_20_30 = PSI_tab.index[(np.abs(0.5 - PSI_tab.mean(axis = 1))>0.2) & (np.abs(0.5 - PSI_tab.mean(axis = 1))<=0.30)] 
    exons_30_40 = PSI_tab.index[(np.abs(0.5 - PSI_tab.mean(axis = 1))>0.1) & (np.abs(0.5 - PSI_tab.mean(axis = 1))<=0.20)] 
    exons_40_50 = PSI_tab.index[(np.abs(0.5 - PSI_tab.mean(axis = 1))<=0.1)] 

    exons_obs_50_60 = PSI_tab.index[(PSI_tab.isna().mean(axis=1) <= 0.5) & (PSI_tab.isna().mean(axis=1) > 0.4)]
    exons_obs_60_70 = PSI_tab.index[(PSI_tab.isna().mean(axis=1) <= 0.4) & (PSI_tab.isna().mean(axis=1) > 0.3)]
    exons_obs_70_80 = PSI_tab.index[(PSI_tab.isna().mean(axis=1) <= 0.3) & (PSI_tab.isna().mean(axis=1) > 0.2)]
    exons_obs_80_90 = PSI_tab.index[(PSI_tab.isna().mean(axis=1) <= 0.2) & (PSI_tab.isna().mean(axis=1) > 0.1)]
    exons_obs_90_100 = PSI_tab.index[(PSI_tab.isna().mean(axis=1) <= 0.1)]

    list1 = [exons_05_10, exons_10_20, exons_20_30, exons_30_40, exons_40_50]
    list2 = [exons_obs_50_60, exons_obs_60_70, exons_obs_70_80, exons_obs_80_90, exons_obs_90_100]


    exon_out_list = []
    C_score_list = []

    for lista_1 in list1:
        for lista_2 in list2:
            combination = lista_1 & lista_2
            if len(combination) > 0:
                exon_out, C_scores = make_mock_C_scores(norm_PSI, Ws, lista_1&lista_2, total_cells, mock=mock)
                exon_out_list.append(exon_out)
                C_score_list.append(C_scores)


    psi_key = ['psi_05_10', 'psi_10_20', 'psi_20_30', 'psi_30_40', 'psi_40_50']
    obs_key = ['obs_50_60', 'obs_60_70', 'obs_70_80', 'obs_80_90', 'obs_90_100']

    counter = 0
    mock_dict = {}
    for pk in psi_key:
        obs_dict = {}
        for ok in obs_key:

            #fit_alpha, fit_loc, fit_beta = st.gamma.fit(C_score_list[counter])
            #random_data = st.gamma.rvs(fit_alpha, loc=fit_loc, scale=fit_beta, size=1000000)
            random_data = C_score_list[counter]
            obs_dict.update({ok:random_data})
            counter += 1
        mock_dict.update({pk:obs_dict})

    return mock_dict



#######################


def get_C_score_pval_gamma(PSI_tab, norm_PSI, Ws, exon_list, total_cells, mock_dict):
    
    
    exon_out_list = []
    C_list = []
    p_list = []
    
    for exon in tqdm(exon_list):
        psi_mean = PSI_tab.loc[exon].mean()
        obs_mean = PSI_tab.loc[exon].isna().mean()
        
        
        
        exon_df = norm_PSI.loc[exon] # to make things faster
        exon_score = get_C(exon_df, Ws)
        if exon_score >= 0:
            C_list.append(exon_score)
            exon_out_list.append(exon)
            
            if (np.abs(0.5 - psi_mean) > 0.4) and (np.abs(0.5 - psi_mean) <= 0.45):
                pk = 'psi_05_10'
            elif (np.abs(0.5 - psi_mean) > 0.3) and (np.abs(0.5 - psi_mean) <= 0.4):
                pk = 'psi_10_20'
            elif (np.abs(0.5 - psi_mean) > 0.2) and (np.abs(0.5 - psi_mean) <= 0.3):
                pk = 'psi_20_30'
            elif (np.abs(0.5 - psi_mean) > 0.1) and (np.abs(0.5 - psi_mean) <= 0.2):
                pk = 'psi_30_40'
            elif (np.abs(0.5 - psi_mean) <= 0.1):
                pk = 'psi_40_50'
                
            if (obs_mean <= 0.5) and (obs_mean > 0.4):
                ok = 'obs_50_60'
            elif (obs_mean <= 0.4) and (obs_mean > 0.3):
                ok = 'obs_60_70'
            elif (obs_mean <= 0.3) and (obs_mean > 0.2):
                ok = 'obs_70_80'
            elif (obs_mean <= 0.2) and (obs_mean > 0.1):
                ok = 'obs_80_90'
            elif (obs_mean <= 0.1):
                ok = 'obs_90_100'
                
            
            random_data = mock_dict[pk][ok]
            
            x = np.sum(random_data > exon_score)
            n = len(random_data)
            pv = (x+1)/(n+1)
            
            p_list.append(pv)
            
    pval_df = pd.DataFrame()
    pval_df['C_score'] = C_list
    pval_df['pval'] = p_list
    pval_df.index = exon_out_list
    return pval_df

     
######################

print('')
print('#######################')
print('Working on Chen')

W_chen = get_distance_matrix(chen_pca_clust, k=244)
chen_norm_PSI = get_signature_matrix(chen_PSI)

int_genes, int_exons = spu.get_int_events(chen_PSI, chen_mrna_counts, 0.05)
observed_exons_1 = chen_PSI.index[chen_PSI[chen_pca_clust.index[chen_pca_clust.AC==0]].isna().mean(axis=1) <= (1-0.5)]
observed_exons_2 = chen_PSI.index[chen_PSI[chen_pca_clust.index[chen_pca_clust.AC==1]].isna().mean(axis=1) <= (1-0.5)]
observed_exons_3 = chen_PSI.index[chen_PSI[chen_pca_clust.index[chen_pca_clust.AC==2]].isna().mean(axis=1) <= (1-0.5)]
observed_exons_4 = chen_PSI.index[chen_PSI[chen_pca_clust.index[chen_pca_clust.AC==3]].isna().mean(axis=1) <= (1-0.5)]
observed_exons_5 = chen_PSI.index[chen_PSI[chen_pca_clust.index[chen_pca_clust.AC==4]].isna().mean(axis=1) <= (1-0.5)]

test_exons = []
for exon in int_exons:
    exon_counts = 0
    exon_counts += (exon in observed_exons_1)
    exon_counts += (exon in observed_exons_2)
    exon_counts += (exon in observed_exons_3)
    exon_counts += (exon in observed_exons_4)
    exon_counts += (exon in observed_exons_5)
    if exon_counts >= 3:
        test_exons.append(exon)



chen_mock_dict = get_mock_dict(chen_PSI.loc[test_exons], chen_norm_PSI.loc[test_exons], W_chen, mock=20000)

pgamma_chen = get_C_score_pval_gamma(chen_PSI.loc[test_exons], chen_norm_PSI.loc[test_exons], 
                                W_chen, test_exons, 488, chen_mock_dict)

pgamma_chen.to_csv('autocorrelation_results/chen.autocorrelation.tab', sep='\t', header=True, index=True)



#############################

print('')
print('#######################')
print('Working on Trapnell')


trapnell_norm_PSI = get_signature_matrix(trapnell_PSI)
W_trapnell = get_distance_matrix(trapnell_pca_clust, k=115)

int_genes, int_exons = spu.get_int_events(trapnell_PSI, trapnell_mrna_counts, 0.05)
observed_exons_1 = trapnell_PSI.index[trapnell_PSI[trapnell_pca_clust.index[trapnell_pca_clust.AC==0]].isna().mean(axis=1) <= (1-0.5)]
observed_exons_2 = trapnell_PSI.index[trapnell_PSI[trapnell_pca_clust.index[trapnell_pca_clust.AC==1]].isna().mean(axis=1) <= (1-0.5)]
observed_exons_3 = trapnell_PSI.index[trapnell_PSI[trapnell_pca_clust.index[trapnell_pca_clust.AC==2]].isna().mean(axis=1) <= (1-0.5)]
observed_exons_4 = trapnell_PSI.index[trapnell_PSI[trapnell_pca_clust.index[trapnell_pca_clust.AC==3]].isna().mean(axis=1) <= (1-0.5)]

test_exons = []
for exon in int_exons:
    exon_counts = 0
    exon_counts += (exon in observed_exons_1)
    exon_counts += (exon in observed_exons_2)
    exon_counts += (exon in observed_exons_3)
    exon_counts += (exon in observed_exons_4)
    if exon_counts >= 3:
        test_exons.append(exon)

trapnell_mock_dict = get_mock_dict(trapnell_PSI.loc[test_exons], trapnell_norm_PSI.loc[test_exons], W_trapnell, mock=20000)

pgamma_trapnell = get_C_score_pval_gamma(trapnell_PSI.loc[test_exons], trapnell_norm_PSI.loc[test_exons], 
                                W_trapnell, test_exons, len(trapnell_PSI.columns), trapnell_mock_dict)

pgamma_trapnell.to_csv('autocorrelation_results/trapnell.autocorrelation.tab', sep='\t', header=True, index=True)





########################

print('')
print('#######################')
print('Working on Song')

song_norm_PSI = get_signature_matrix(song_PSI)
W_song = get_distance_matrix(song_pca_clust, k=101)

int_genes, int_exons = spu.get_int_events(song_PSI, song_mrna_counts, 0.05)
observed_exons_1 = song_PSI.index[song_PSI[song_pca_clust.index[song_pca_clust.AC==0]].isna().mean(axis=1) <= (1-0.5)]
observed_exons_2 = song_PSI.index[song_PSI[song_pca_clust.index[song_pca_clust.AC==1]].isna().mean(axis=1) <= (1-0.5)]
observed_exons_3 = song_PSI.index[song_PSI[song_pca_clust.index[song_pca_clust.AC==2]].isna().mean(axis=1) <= (1-0.5)]
observed_exons = observed_exons_1 & observed_exons_2 & observed_exons_3
test_exons = [x for x in observed_exons if x in int_exons]

song_mock_dict = get_mock_dict(song_PSI.loc[test_exons], song_norm_PSI.loc[test_exons], W_song, mock=20000)

pgamma_song = get_C_score_pval_gamma(song_PSI.loc[test_exons], song_norm_PSI.loc[test_exons], 
                                W_song, test_exons, len(song_PSI.columns), song_mock_dict)

pgamma_song.to_csv('autocorrelation_results/song.autocorrelation.tab', sep='\t', header=True, index=True)


########################

















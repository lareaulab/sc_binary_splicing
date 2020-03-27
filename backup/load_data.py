import numpy as np
import pandas as pd
import os
import matplotlib.cm as cm
from matplotlib import pyplot as plt
from scipy import stats as st
import seaborn as sns

from IPython.core.pylabtools import figsize

import numpy.random as r
from pylab import *
from matplotlib.gridspec import GridSpec

import sys
sys.path.insert(0, '../../utils')
import splicing_utils as spu
from splicing_utils import *
import single_cell_plots as scp
from single_cell_plots import *

plt.rcParams["axes.edgecolor"] = "black"
plt.rcParams["axes.linewidth"] = 1
plt.rcParams["axes.facecolor"] = 'white'

import matplotlib as mpl
import numpy as np
from matplotlib import pyplot as plt

from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering

mpl.rcParams["mathtext.fontset"] = "stix"


# PSI tables
chen_PSI = pd.read_csv('../../../data/chen/processed_tables/chen.skipped_exons_psi.tab', sep='\t', index_col=0)
song_PSI = pd.read_csv('../../../data/song/processed_tables/song.skipped_exons_psi.tab', sep='\t',  index_col=0)
das_PSI = pd.read_csv('../../../data/das/processed_tables/das.skipped_exons_psi.tab', sep='\t',  index_col=0)
trapnell_PSI = pd.read_csv('../../../data/trapnell/processed_tables/trapnell.skipped_exons_psi.tab', 
                    sep='\t',  index_col=0)
lescroart_PSI = pd.read_csv('../../../data/lescroart/processed_tables/lescroart.skipped_exons_psi.tab', 
                     sep='\t',  index_col=0)
shalek_PSI = pd.read_csv('../../../data/shalek/processed_tables/shalek.skipped_exons_psi.tab', sep='\t', index_col=0)

# SJ read tables
chen_read_counts = pd.read_csv('../../../data/chen/processed_tables/chen.skipped_exons_SJreads.tab', sep='\t', index_col=0)
song_read_counts = pd.read_csv('../../../data/song/processed_tables/song.skipped_exons_SJreads.tab', sep='\t', index_col=0)
das_read_counts = pd.read_csv('../../../data/das/processed_tables/das.skipped_exons_SJreads.tab', sep='\t',  index_col=0)
trapnell_read_counts = pd.read_csv('../../../data/trapnell/processed_tables/trapnell.skipped_exons_SJreads.tab', 
                            sep='\t', index_col=0)
lescroart_read_counts = pd.read_csv('../../../data/lescroart/processed_tables/lescroart.skipped_exons_SJreads.tab', 
                             sep='\t', index_col=0)
shalek_read_counts = pd.read_csv('../../../data/shalek/processed_tables/shalek.skipped_exons_SJreads.tab', 
                          sep='\t', index_col=0)

# TPM tables
chen_tpm_tab = pd.read_csv('../../../data/chen/chen.tpm.gene_symbols.tab', sep='\t', index_col=0)
song_tpm_tab = pd.read_csv('../../../data/song/song.tpm.gene_symbols.tab', sep='\t', index_col=0)
trapnell_tpm_tab = pd.read_csv('../../../data/trapnell/trapnell.tpm.gene_symbols.tab', sep='\t', index_col=0)
lescroart_tpm_tab = pd.read_csv('../../../data/lescroart/lescroart.tpm.gene_symbols.tab', sep='\t', index_col=0)
shalek_tpm_tab = pd.read_csv('../../../data/shalek/shalek.tpm.gene_symbols.tab', sep='\t', index_col=0)
das_tpm_tab = pd.read_csv('../../../data/das/das.tpm.gene_symbols.tab', sep='\t', index_col=0)



chen_counts_tab = pd.read_csv('../../../data/chen/processed_tables/chen.rsemCounts.gene_symbols.tab', 
                       sep='\t', index_col=0)
song_counts_tab = pd.read_csv('../../../data/song/processed_tables/song.rsemCounts.gene_symbols.tab', 
                       sep='\t', index_col=0)
das_counts_tab = pd.read_csv('../../../data/das/processed_tables/das.rsemCounts.gene_symbols.tab', 
                      sep='\t', index_col=0)
trapnell_counts_tab = pd.read_csv('../../../data/trapnell/processed_tables/trapnell.rsemCounts.gene_symbols.tab', 
                           sep='\t', index_col=0)
lescroart_counts_tab = pd.read_csv('../../../data/lescroart/processed_tables/lescroart.rsemCounts.gene_symbols.tab', 
                            sep='\t', index_col=0)
shalek_counts_tab = pd.read_csv('../../../data/shalek/processed_tables/shalek.rsemCounts.gene_symbols.tab', 
                         sep='\t', index_col=0)

# mRNA tables
chen_mrna_counts = pd.read_csv('../../../data/chen/processed_tables/chen.mrna_counts.tab', sep='\t', index_col=0)
song_mrna_counts = pd.read_csv('../../../data/song/processed_tables/song.mrna_counts.tab', sep='\t', index_col=0)
das_mrna_counts = pd.read_csv('../../../data/das/processed_tables/das.mrna_counts.tab', sep='\t', index_col=0)
trapnell_mrna_counts = pd.read_csv('../../../data/trapnell/processed_tables/trapnell.mrna_counts.tab', 
                            sep='\t', index_col=0)
lescroart_mrna_counts = pd.read_csv('../../../data/lescroart/processed_tables/lescroart.mrna_counts.tab', 
                             sep='\t', index_col=0)
shalek_mrna_counts = pd.read_csv('../../../data/shalek/processed_tables/shalek.mrna_counts.tab', 
                          sep='\t', index_col=0)

# mRNA per event rables
mrna_per_event_chen = pd.read_csv('../../../data/chen/processed_tables/chen.mrna_counts_per_event.tab', sep='\t', index_col=0)
mrna_per_event_song = pd.read_csv('../../../data/song/processed_tables/song.mrna_counts_per_event.tab', sep='\t', index_col=0)
mrna_per_event_das = pd.read_csv('../../../data/das/processed_tables/das.mrna_counts_per_event.tab', sep='\t', index_col=0)
mrna_per_event_trapnell = pd.read_csv('../../../data/trapnell/processed_tables/trapnell.mrna_counts_per_event.tab', 
                            sep='\t', index_col=0)
mrna_per_event_lescroart = pd.read_csv('../../../data/lescroart/processed_tables/lescroart.mrna_counts_per_event.tab', 
                             sep='\t', index_col=0)
mrna_per_event_shalek = pd.read_csv('../../../data/shalek/processed_tables/shalek.mrna_counts_per_event.tab', 
                          sep='\t', index_col=0)

# read coverage tables
chen_coverage_tab = pd.read_csv('../../../data/chen/processed_tables/chen.read_coverage.tab', 
                          sep='\t', index_col=0)
song_coverage_tab = pd.read_csv('../../../data/song/processed_tables/song.read_coverage.tab', 
                          sep='\t', index_col=0)
trapnell_coverage_tab = pd.read_csv('../../../data/trapnell/processed_tables/trapnell.read_coverage.tab', 
                          sep='\t', index_col=0)
lescroart_coverage_tab = pd.read_csv('../../../data/lescroart/processed_tables/lescroart.read_coverage.tab', 
                          sep='\t', index_col=0)
das_coverage_tab = pd.read_csv('../../../data/das/processed_tables/das.read_coverage.tab', 
                          sep='\t', index_col=0)
shalek_coverage_tab = pd.read_csv('../../../data/shalek/processed_tables/shalek.read_coverage.tab', 
                          sep='\t', index_col=0)


chen_pca = pd.read_csv('../../../data/chen/chen.pca.tab', sep='\t', index_col=0)
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

song_pca = pd.read_csv('../../../data/song/song.pca_top_fano.tab', sep='\t', index_col=0)
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

trapnell_pca = pd.read_csv('../../../data/trapnell/trapnell.pca.tab', sep='\t', index_col=0)
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
        
lescroart_pca = pd.read_csv('../../../data/lescroart/lescroart.pca_meta.tab', sep='\t', index_col = 0)
lescroart_index = [x for x in lescroart_pca.index if x in mrna_per_event_lescroart.columns]

lescroart_pca = lescroart_pca.loc[lescroart_index]
lescroart_PSI = lescroart_PSI[lescroart_index]
mrna_per_event_lescroart = mrna_per_event_lescroart[lescroart_index]
lescroart_read_counts = lescroart_read_counts[lescroart_index]
lescroart_coverage_tab = lescroart_coverage_tab.loc[lescroart_index]

lescroart_E6 = lescroart_pca.loc[lescroart_pca.cell_type=='E6.75'].index
lescroart_E7 = lescroart_pca.loc[lescroart_pca.cell_type=='E7.25'].index

def process_subpop(subpop, psi, mrna, mrna_per_event, reads, cj, psi_min = 0.2, mrna_min=10, reads_min = 0, cell_min = 0.5, nbins=11):

    int_genes, int_exons = spu.get_int_events(psi[subpop], mrna[subpop], psi_min)
    #print(len(int_genes))
    int_exons = [x for x in int_exons if x in mrna_per_event.index]
    PSI_filtered, PSI_mrna_filtered, good_exons, mrna_filtered, reads_filtered = filter_psi(psi[subpop], int_exons, 
                                                                     mrna_per_event[subpop], cj.loc[subpop], 
                                                                     reads[subpop], mrna_min, reads_min = reads_min,
                                                                                            cell_min=cell_min)


    good_cells = PSI_filtered.dropna(axis=1, how='all').columns
    good_subpop = [x for x in subpop if x in good_cells]
    PSI_good = PSI_filtered[good_cells]

    hist_complete, hist_intermediate = scp.get_bins_table2(PSI_filtered[good_subpop], mrna_filtered[good_subpop], nbins)
    hist_complete_exp, hist_intermediate_exp = scp.get_bins_table(PSI_filtered[good_subpop], mrna_filtered[good_subpop])

    
    return PSI_filtered, good_exons, mrna_filtered, reads_filtered, hist_complete, hist_complete_exp
    


##############################

###################
# Chen
ac = AgglomerativeClustering(n_clusters=5)
ac_clusters = ac.fit_predict(chen_pca[['PC1', 'PC2']])

# figsize(6,4)
# plt.scatter(chen_pca.PC1, chen_pca.PC2, c=ac_clusters)
# plt.show()

chen_pca_clust = chen_pca.copy()
chen_pca_clust['AC'] = ac_clusters

chen_clust_filter = []
for cluster in chen_pca_clust.groupby('AC')['pseudotime'].mean().sort_values().index:
    clust_subpop = chen_pca_clust.index[chen_pca_clust.AC == cluster]
    
    chen_filter = process_subpop(clust_subpop, chen_PSI, chen_mrna_counts, mrna_per_event_chen, 
                                 chen_read_counts, chen_coverage_tab['SJ_coverage'], 0.1, 10, 0, cell_min=0.5)
    
    chen_clust_filter.append(chen_filter)
    
###################
# Song
ac = AgglomerativeClustering(n_clusters=3)
ac_clusters = ac.fit_predict(song_pca[['PC1', 'PC2']])

# figsize(6,4)
# plt.scatter(song_pca.PC1, song_pca.PC2, c=ac_clusters)
# plt.show()

song_pca_clust = song_pca.copy()
song_pca_clust['AC'] = ac_clusters

song_clust_filter = []
for cluster in song_pca_clust.groupby('AC')['pseudotime'].mean().sort_values().index:
    clust_subpop = song_pca_clust.index[song_pca_clust.AC == cluster]
    
    song_filter = process_subpop(clust_subpop, song_PSI, song_mrna_counts, mrna_per_event_song, 
                                 song_read_counts, song_coverage_tab['SJ_coverage'], 0.1, 10, 0, cell_min=0.5)
    
    song_clust_filter.append(song_filter)
    
    
###################
# Lescroart
ac = AgglomerativeClustering(n_clusters=3)
ac_clusters = ac.fit_predict(lescroart_pca[['PC1', 'PC2']])

# figsize(6,4)
# plt.scatter(lescroart_pca.PC1, lescroart_pca.PC2, c=ac_clusters)
# plt.show()

lescroart_pca_clust = lescroart_pca.copy()
lescroart_pca_clust['AC'] = ac_clusters


lescroart_clust_filter = []
for cluster in range(3):
    clust_subpop = lescroart_pca_clust.index[lescroart_pca_clust.AC == cluster]
    
    lescroart_filter = process_subpop(clust_subpop, lescroart_PSI, lescroart_mrna_counts, mrna_per_event_lescroart, 
                                 lescroart_read_counts, lescroart_coverage_tab['SJ_coverage'], 0.1, 10, 0, cell_min=0.5)
    
    lescroart_clust_filter.append(lescroart_filter)
    
    
###################
# Trapnell
ac = AgglomerativeClustering(n_clusters=4)
ac_clusters = ac.fit_predict(trapnell_pca[['PC1', 'PC2']])

# figsize(6,4)
# plt.scatter(trapnell_pca.PC1, trapnell_pca.PC2, c=ac_clusters)
# plt.show()

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
import rpy2
import rpy2.robjects.packages as rpackages
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri as rpyn
from statsmodels.stats.multitest import multipletests
dt = rpy2.robjects.packages.importr('diptest')

from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering

from scipy.special import logit
from scipy.special import expit

from sklearn.metrics import adjusted_rand_score

from scipy.stats import chisquare
from scipy.stats import hypergeom
from scipy.stats import f_oneway
from scipy.stats import kruskal

from scipy.stats import combine_pvalues
from scipy.stats import probplot
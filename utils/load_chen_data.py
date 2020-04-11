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
sys.path.insert(0, '.')
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
mpl.rcParams['pdf.fonttype'] = 42

data_dir = '/mnt/c/Users/ferna/Desktop/SingleCell/data/'

# PSI tables
chen_PSI = pd.read_csv(data_dir + 'chen/processed_tables/chen.skipped_exons_psi.tab', sep='\t', index_col=0)

# SJ read tables
chen_read_counts = pd.read_csv(data_dir + 'chen/processed_tables/chen.skipped_exons_SJreads.tab', sep='\t', index_col=0)

# TPM tables
# chen_tpm_tab = pd.read_csv(data_dir + 'chen/chen.tpm.gene_symbols.tab', sep='\t', index_col=0)

#chen_counts_tab = pd.read_csv(data_dir + 'chen/processed_tables/chen.rsemCounts.gene_symbols.tab', 
#                       sep='\t', index_col=0)
# mRNA tables
chen_mrna_counts = pd.read_csv(data_dir + 'chen/processed_tables/chen.mrna_counts.tab', sep='\t', index_col=0)

# mRNA per event rables
mrna_per_event_chen = pd.read_csv(data_dir + 'chen/processed_tables/chen.mrna_counts_per_event.tab', sep='\t', index_col=0)

# read coverage tables
chen_coverage_tab = pd.read_csv(data_dir + 'chen/processed_tables/chen.read_coverage.tab', 
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
    

#!/usr/bin/env python
# coding: utf-8

# # Simulations 2
# 
# In these simulations, we test the effect of different capture efficiencies in the observed binarity of the data. We simulate single cell splicing with a fixed $\Psi = 0.5$ for every exon, and then simulate mRNA capture and sequencing with different captur efficiency levels $c\in\{0.01, 0.012, ..., 0.1\}$. We repeated the simulations 50 times for each capture efficiency value.
# 
# We then estimated the observed $\Psi$ for each exon in each cell. We ranked the exon by the number of reads that cover them, and estimated the fraction of binary ($\Psi = 0/1$) observations. We also estimated the deviation of the estimated $\Psi$ from the real $\Psi = 0.5$.

# In[1]:


import numpy as np
import pandas as pd
import os
from matplotlib import pyplot as plt
from numpy import random as r

from pylab import *
from IPython.core.pylabtools import figsize

import seaborn as sns
import matplotlib.cm as cm

import sys
sys.path.insert(0, '/mnt/lareaulab/cfbuenabadn/SingleCell/utils')
import splicing_utils as spu
import single_cell_plots as scp
from single_cell_plots import *


# In[2]:


def process_symsim(tabla, gene_len, n=300, g=500, read_len=50):
    
    gene_len_i = gene_len[:g]
    gene_len_e = gene_len[g:]
    
    biases_i = (4*(read_len-1))/gene_len_i
    biases_e = (2*(read_len-1))/gene_len_e
    
    i = tabla.iloc[:g]
    e = tabla.iloc[g:].reset_index(drop=True)
    
    cuentas_i = pd.DataFrame(np.array([r.binomial(i.loc[x], biases_i[x]) for x in i.index]))
    cuentas_e = pd.DataFrame(np.array([r.binomial(e.loc[x], biases_e[x]) for x in e.index]))
    
    counts = cuentas_i + 2*cuentas_e
    psi = cuentas_i/counts
    
    return counts, psi#, i, e




def get_bimodal(obs_table, gene_len):
    counts_table, psi_table = process_symsim(obs_table, gene_len)
    counts_table.columns = obs_table.columns
    psi_table.columns = obs_table.columns

    counts_table.index = ['ase_'+str(x) for x in range(1, 501)]
    psi_table.index = ['ase_'+str(x) for x in range(1, 501)]
    
    sorted_mean_expression = counts_table.mean(axis=1).sort_values()
    sorted_ase = sorted_mean_expression.index
    
    psi_bimodal = (psi_table.loc[sorted_ase] <= 0.01) | (psi_table.loc[sorted_ase] >= 0.99)
    
    sorted_error = np.sum(psi_bimodal, axis=1) / ((~psi_table.loc[sorted_ase].isnull()).sum(axis=1))
    
    return sorted_mean_expression, sorted_error


def get_error(obs_table, gene_len):
    counts_table, psi_table = process_symsim(obs_table, gene_len)
    counts_table.columns = obs_table.columns
    psi_table.columns = obs_table.columns

    counts_table.index = ['ase_'+str(x) for x in range(1, 501)]
    psi_table.index = ['ase_'+str(x) for x in range(1, 501)]
    
    sorted_mean_expression = counts_table.mean(axis=1).sort_values()
    sorted_ase = sorted_mean_expression.index
    sorted_error = np.sum(np.abs(psi_table.loc[sorted_ase] - 0.5), axis=1) / ((~psi_table.loc[sorted_ase].isnull()).sum(axis=1))
    
    return sorted_mean_expression, sorted_error


# In[12]:


sim_dir = '/mnt/lareaulab/cfbuenabadn/SingleCell/Results/simulations/'

hm = pd.DataFrame()
out_hm = 0

for sim in range(1,51):
    sim_s = str(sim)
    gene_len_str = sim_dir + 'sim_tables11/gene_length.sim_' + sim_s + '.txt'
    gene_len = np.array([int(x.rstrip()) for x in open(gene_len_str).readlines()])

    for cap in np.arange(0.01, 0.101, step=0.002):

        cap_s = str(round(cap, 3))

        suffix = '.sim_' + sim_s + '_cap_' + cap_s

        table = pd.read_csv(sim_dir + 'sim_tables11/observed_counts' + suffix + '.tab', sep='\t', 
                          names = ['cell_'+str(x) for x in range(1, 301)])

        sorted_mean_expression, sorted_error = get_bimodal(table, gene_len)

        hm['cap_' + cap_s] = list(sorted_error)[::-1]
    out_hm += hm
    print(sim)
    
out_hm = out_hm/50


# In[13]:


out_hm.columns = [x.split('_')[1] for x in hm.columns]
out_hm.index = range(1,501)

out_hm.to_csv('simulations_2_binary.tab', sep='\t', index=True, header=True)


# In[21]:


#fig = plt.figure(figsize=(12,8))

#figure_title = "Capture efficiency and binary values"
#ax  = plt.subplot(1,1,1)

#plt.text(0.5, 1.03, figure_title,
#         horizontalalignment='center',
#         fontsize=28,
#         transform = ax.transAxes)
#sns.heatmap(out_hm.dropna(), cmap='bwr')
#
#plt.xlabel('Capture efficiency', fontsize=24)
#plt.ylabel('ASE ranked by SJ', fontsize=24)

#plt.show()


# In[3]:


sim_dir = '/mnt/lareaulab/cfbuenabadn/SingleCell/Results/simulations/'

hm = pd.DataFrame()
error_hm = 0

for sim in range(1,51):
    sim_s = str(sim)
    gene_len_str = sim_dir + 'sim_tables11/gene_length.sim_' + sim_s + '.txt'
    gene_len = np.array([int(x.rstrip()) for x in open(gene_len_str).readlines()])

    for cap in np.arange(0.01, 0.101, step=0.002):

        cap_s = str(round(cap, 3))

        suffix = '.sim_' + sim_s + '_cap_' + cap_s

        table = pd.read_csv(sim_dir + 'sim_tables11/observed_counts' + suffix + '.tab', sep='\t', 
                          names = ['cell_'+str(x) for x in range(1, 301)])

        sorted_mean_expression, sorted_error = get_error(table, gene_len)

        hm['cap_' + cap_s] = list(sorted_error)[::-1]
    error_hm += hm
    print(sim)
    
error_hm = error_hm/50


# In[ ]:


error_hm.columns = [x.split('_')[1] for x in hm.columns]
error_hm.index = range(1,501)

error_hm.to_csv('simulations_2_error.tab', sep='\t', index=True, header=True)


# In[ ]:





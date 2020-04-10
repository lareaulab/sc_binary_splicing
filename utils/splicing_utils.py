import numpy as np
import pandas as pd
import subprocess as sp
import os
from scipy.stats import gaussian_kde
from numpy import random as r
import matplotlib.pyplot as plt


def get_psi_table(SJ_table_name, minJR=5, minCell=20, drop_duplicates = False):
    '''
    
    This functions splits this table into one individual for
    each junction type. It additionally creates a PSI table
    based on PSI = (I1 + I2) / ((I1 + I2) + 2*SE)
    
    Input:
    - SJ_counts_table is a pandas dataframe, with rows corresponding
      to individual splicing junctions, and columns corresponding to
      single cells.
    
      The format of the index is:
      Gene_X_SJ
      Where Gene correspond to the gene where the splicing event
      is present. X is a number assigned to one particulat ASE. 
      NMD events have additionally _nmdSE_ between Gene and X.
      SJ corresponds to I1 (included 1), I2 (included 2) or SE
      (skipped exon).
    
    - minCell (int) is the minimum number of cells that are required
      to have at least minJR reads.
    
    - minJR (int) determines the minimum number of reads required on
      at least minCell cells for a junction to be acceptable.
    
    Output:
    - I1_counts (DF) Counts for one Included SJ
    - I2_counts (DF) Counts for the other Included SJ 
    - SE_counts (DF) Counts for Skipped Exon SJ 
    - PSI_table (DF) As in name.
    - total_counts (DF) SE + (I1 + I2)/2
    
    '''
    
    if drop_duplicates:
        SJ_counts_table = pd.read_csv(SJ_table_name, sep='\t', index_col=0).drop_duplicates('last')
    else:
        SJ_counts_table = pd.read_csv(SJ_table_name, sep='\t', index_col=0)
    
    # Select dataframes for each junction type
    I1 = [x for x in SJ_counts_table.index if '_I1' in x]
    I2 = [x for x in SJ_counts_table.index if '_I2' in x]
    SE = [x for x in SJ_counts_table.index if '_SE' in x]

    # Filter out junctions with very low coverage
    I1_filt = (SJ_counts_table.loc[I1] > minJR).sum(axis=1) > minCell
    I2_filt = (SJ_counts_table.loc[I2] > minJR).sum(axis=1) > minCell
    SE_filt = (SJ_counts_table.loc[SE] > minJR).sum(axis=1) > minCell

    # Remove junction extension from index (to be able to merge them)
    I1_idx = [x[:-3] for x in I1_filt.index if I1_filt[x]]
    I2_idx = [x[:-3] for x in I2_filt.index if I2_filt[x]]
    SE_idx = [x[:-3] for x in SE_filt.index if SE_filt[x]]

    # ASE events that pass the minimum requirements in the three SJ
    filtered_index = [x for x in I1_idx if ((x in SE_idx) and (x in I2_idx))]
    I1_counts = SJ_counts_table.loc[[x + '_I1' for x in filtered_index]]
    I2_counts = SJ_counts_table.loc[[x + '_I2' for x in filtered_index]]
    SE_counts = SJ_counts_table.loc[[x + '_SE' for x in filtered_index]]

    # Correct index
    I1_counts.index = filtered_index
    I2_counts.index = filtered_index
    SE_counts.index = filtered_index

    # Create PSI table and total counts table
    PSI_table = (I1_counts + I2_counts) /(2*SE_counts + I1_counts + I2_counts)
    PSI_table_I = (I1_counts) /(SE_counts + I1_counts)
    total_counts = SE_counts + I1_counts + I2_counts
    total_counts_I = SE_counts + I1_counts
    
    
    clean_SJ2 = np.intersect1d(np.array(total_counts.drop_duplicates().index), np.array(PSI_table.drop_duplicates().index))
    clean_SJ1 = np.intersect1d(np.array(total_counts_I.drop_duplicates().index), np.array(PSI_table_I.drop_duplicates().index))
    
    PSI_table = PSI_table.loc[clean_SJ2]
    total_counts = total_counts.loc[clean_SJ2]
    
    PSI_table_I = PSI_table_I.loc[clean_SJ1]
    total_counts_I = total_counts_I.loc[clean_SJ1]
    
    return I1_counts, I2_counts, SE_counts, PSI_table, total_counts, PSI_table_I, total_counts_I


def process_symsim_true(table_name, n=500, g=1500):
    
    '''
    input:
    * table_name (str): path to the table with isoform counts. 
      Should have a (2g,n) dimension.
    * n (int): number of cells in the simulation (columns)
    * g (int): number of genes in the simulation. Isoforms are 2g.
    
    output:
    * counts (df): dataframe (g,n) with the total mRNA or read counts 
      for each gene (sum of included and excluded counts).
    * psi (df): dataframe (g,n) with ratio of included / total counts. NA
      for cases where total counts == 0.
      
    Description:    
      Get the ratio of isoform read or mRNA counts from SymSim output.
      No splicing junction downsampling.
    
      Assumes that the original table's shape is (n, 2g). The first 
      g rows correspond to the "included" isoforms. The latter g rows
      correspond to the "excluded" isoform (order is preserved).
    
      Originally used for reads before I realized the bias SymSim introduces 
      when gene length is set as constant.
    
      Useful for getting the true PSI of ASE.
    '''
    
    # read table
    tabla = pd.read_csv(table_name, sep='\t', names = ['cell_' + str(i) for i in range(1, n+1)])
    
    # select counts for each isoform
    i = tabla.iloc[:g]
    e = tabla.iloc[g:].reset_index(drop=True)
    
    # calculate total counts and PSI
    counts = i + e
    psi = i/counts
    
    return counts, psi

def process_symsim_observed(table_name, gene_len, n=500, g=1500, read_len=50):
    
    gene_len_i = gene_len[:g]
    gene_len_e = gene_len[g:]
    
    biases_i = (4*(read_len-1))/gene_len_i
    biases_e = (2*(read_len-1))/gene_len_e
    
    tabla = pd.read_csv(table_name, sep='\t', names = ['cell_' + str(i) for i in range(1, n+1)])
    
    i = tabla.iloc[:g]
    e = tabla.iloc[g:].reset_index(drop=True)
    
    cuentas_i = pd.DataFrame(np.array([r.binomial(i.loc[x], biases_i[x]) for x in i.index]))
    cuentas_e = pd.DataFrame(np.array([r.binomial(e.loc[x], biases_e[x]) for x in e.index]))

    counts = cuentas_i + 2*cuentas_e
    psi = cuentas_i/counts
    
    return counts, psi#, i, e


def best_cells(counts_df, psi_df, fraction=2):
    half_cells = int(np.floor(len(counts_df.columns)/fraction))
    filtered_counts = pd.DataFrame()
    filtered_psi = pd.DataFrame()
    for idx in counts_df.index:
        top_cells = counts_df.loc[idx].sort_values().index[-half_cells:]
        filtered_counts[idx] = list(counts_df.loc[idx, top_cells])
        filtered_psi[idx] = list(psi_df.loc[idx, top_cells])
        
    return filtered_counts, filtered_psi


def normalize_equation(cell, moda):
    n = np.sum(np.log10(cell) <= moda)
    interval = np.sum(cell.loc[np.log10(cell) <= moda])/np.sum(cell)
    return n/interval
    

def get_transcripts_per_cell(cell, plot_hist, correct_high, bw_method, adjust_high):
    z = gaussian_kde(np.log10(cell), bw_method)
    
    moda = np.arange(-1, 10, 0.1)[z.pdf(np.arange(-1, 10, 0.1)).argmax()]
    
    molecules_in_cell = normalize_equation(cell, moda)
    
    
    if (molecules_in_cell > 1000000) and adjust_high:
        moda = np.arange(1, 10, 0.1)[z.pdf(np.arange(1, 10, 0.1)).argmax()]
    
        molecules_in_cell = normalize_equation(cell, moda)
            
    
    if correct_high and (molecules_in_cell > 1000000):
        
        return 0
        
    if plot_hist:
        plt.figure()
        plt.hist(np.log10(cell), density=True, alpha=0.15)
        plt.plot([moda, moda],
                 [0, z.pdf(moda)[0]], '--', linewidth=3, c='navy')
        plt.plot(np.arange(-1, 5, 0.01), z.pdf(np.arange(-1, 5, 0.01)), c='darkred', linewidth=5)

        plt.show()
    
    return molecules_in_cell
        

def transform_cell(cell, plot_hist, correct_high, bw_method, adjust_high):
    cell_filtered = cell.loc[cell > 0.1]
    molecules_in_cell = get_transcripts_per_cell(cell_filtered, plot_hist, correct_high, bw_method, adjust_high)
    cell_remove_zeros = cell * (cell > 0.1)
    normalize_lysate = molecules_in_cell / 1000000
    cell_transcript_counts = cell_remove_zeros * normalize_lysate
    
    return cell_transcript_counts#, molecules_in_cell


def transform_tpm_to_counts(tpm_dataset, plot_hist = True, correct_high = True, bw_method='scott', adjust_high = False):
    mrna_counts = pd.DataFrame()
    mrna_counts_per_cell = []
    cells = tpm_dataset.columns
    tpm_dataset_filtered = tpm_dataset.loc[tpm_dataset.max(axis=1) > 0.1]
    
    for cell in cells:
        cell_mrna = transform_cell(tpm_dataset_filtered[cell], plot_hist, correct_high, bw_method, adjust_high)
        if all([x == 0 for x in cell_mrna]):
            print('Skipping cell')
            continue
        mrna_counts[cell] = cell_mrna
        
    return mrna_counts


def get_binary_fraction(exon):
    return (np.sum(exon == 0) + np.sum(exon == 1)) / np.sum(~exon.isnull())


def reads_per_mrna(exon_reads, gene_mrna_counts, norm_factor):
    return exon_reads/(norm_factor * gene_mrna_counts)


def get_binary_fraction(exon):
    return (np.sum(exon == 0) + np.sum(exon == 1)) / np.sum(~exon.isnull())

    
def get_rpm_dataset(PSI_table, read_counts, mrna_counts, exon_list):
    
    genes_list = [x.split('_')[0] for x in exon_list]
    cells = PSI_table.columns
    
    cells = [x for x in PSI_table.columns if x in mrna_counts.columns]
    cells = [x for x in cells if x in read_counts.columns]
    
    reads = read_counts.loc[exon_list, cells]
    psi = np.array(PSI_table.loc[exon_list, cells])
    mrna = mrna_counts.loc[genes_list, cells]
    
    reads = np.array(reads) * np.array((mrna > 0).astype(int))
    
    mrna = np.array(mrna)
    
    sr = pd.DataFrame(reads/((1+psi)*mrna))
    
    sr.columns = cells
    sr.index = exon_list

    return sr


# def get_int_events(PSI_table, mrna_counts, min_unimodal = 0.05):
#     int_exons = PSI_table.index[(PSI_table.mean(axis=1) >= min_unimodal) & (PSI_table.mean(axis=1) <= (1- min_unimodal))]
#     int_exons = [x for x in int_exons if x.split('_')[0] in mrna_counts.index]
#     int_genes = sorted(set([x.split('_')[0] for x in int_exons]))
#     return int_genes, int_exons

def get_int_events(PSI_table, mrna_counts, min_unimodal = 0.05):
    int_exons = PSI_table.index[(PSI_table.mean(axis=1) >= min_unimodal) & (PSI_table.mean(axis=1) <= (1- min_unimodal))]
    int_exons = [x for x in int_exons if x.split('_')[0] in mrna_counts.index]
    int_genes = sorted(set([x.split('_')[0] for x in int_exons]))
    return int_genes, int_exons


    
def get_dataset_averages(psi_table, read_counts, mrna_counts, exon_list, stat = 'mean'):
    
    reads_list = []
    mrna_list = []
    binary_list = []
    psi_list = []

    if stat == 'median':
        for e in exon_list:
            reads_list.append(read_counts.loc[e].median())
            mrna_list.append(mrna_counts.loc[e.split('_')[0]].median())
            binary_list.append(get_binary_fraction(psi_table.loc[e]))
            psi_list.append(psi_table.loc[e].median())
        
    else:
        for e in exon_list:
            reads_list.append(read_counts.loc[e].mean())
            mrna_list.append(mrna_counts.loc[e.split('_')[0]].mean())
            binary_list.append(get_binary_fraction(psi_table.loc[e]))
            psi_list.append(psi_table.loc[e].mean())
        
    return reads_list, mrna_list, binary_list, psi_list
    
    
def filter_psi(psi_tab, int_exons, mrna_tab, cj, reads_tab, mrna_min, reads_min=0, cell_min=0.5, filter_cj = True):
    
    '''
    Filter observations and select exons.
    '''
    
    int_exons = [x for x in int_exons if x in mrna_tab.index]
    reads_tab = reads_tab.loc[int_exons, mrna_tab.columns]
    mrna_tab = mrna_tab.loc[int_exons, mrna_tab.columns]
    cj = cj.loc[mrna_tab.columns]
    psi_tab = psi_tab.loc[int_exons, mrna_tab.columns]
    
    if filter_cj:
        cj_factor = 1
    else:
        cj_factor = 0

    cj_filter = (reads_tab.loc[int_exons, cj.index] >= (mrna_min*cj_factor*cj*(1+psi_tab.loc[int_exons, cj.index])))
    #mrna_filter = (mrna_tab.loc[int_exons, cj.index] >= mrna_min) & (reads_tab.loc[int_exons, cj.index] >= reads_min) #&
    ##### reads_min*(1+psi_tab.loc[int_exons, cj.index])
    mrna_filter = (mrna_tab.loc[int_exons, cj.index] >= mrna_min) & (reads_tab.loc[int_exons, cj.index] >= (reads_min*(1+psi_tab.loc[int_exons, cj.index]))) #&
    quality = ((mrna_filter & cj_filter).mean(axis=1) >= cell_min)
    good_exons = quality.loc[quality].index

    filter_tab = (mrna_filter & cj_filter).loc[good_exons]
    psi_tab_mrna = psi_tab.loc[good_exons]

    PSI_filtered = psi_tab_mrna
    PSI_filtered = psi_tab_mrna.mask(~filter_tab)
    PSI_filtered = PSI_filtered.dropna(how='all') # might be unnecessary
    PSI_mrna_filtered = PSI_filtered.loc[PSI_filtered.isnull().mean(axis=1) < (1-cell_min)] # might be unnecessary
    total_exons = len(PSI_mrna_filtered.index)
    
    
    mrna_filtered = mrna_tab.loc[PSI_filtered.index].mask(~filter_tab)
    reads_filtered = reads_tab.loc[PSI_filtered.index].mask(~filter_tab)

    return PSI_filtered, PSI_mrna_filtered, good_exons, mrna_filtered, reads_filtered


def hartigan_test(x_list):
    ht_input = 'c(' + ','.join([str(x) for x in x_list]) + ')'
    ht_test_r = dt.dip_test(x=robjects.r(ht_input))
    ht_test_py = rpyn.rpy2py(ht_test_r)
    ht_statistic = ht_test_py[0][0]
    ht_pval = ht_test_py[1][0]
    return ht_statistic, ht_pval



def filtered_bimodal(subset_list, cell_type, pca_tab, mrna_per_event, 
                     read_counts, quantile = 0.25, psi_lim = 0.25):
    i = 0
    j = 0
    for subset in subset_list:

        for event in subset[0].index:
            
            

            token = subset[0].loc[event].dropna().quantile(quantile) <= psi_lim
            token = token and subset[0].loc[event].dropna().quantile(1-quantile) >= (1 - psi_lim)
            if token:

                plt.hist(subset[0].loc[event].dropna(), range=[0,1], density=True)
                plt.title(event + ', ' + cell_type[j], fontsize=24)
                plt.xlabel(r'$\Psi$', fontsize=24)
                plt.ylabel('frequency', fontsize=24)
                plt.show()

                good_cells = subset[0].loc[event].dropna().index
                all_cells = subset[0].columns
                
                if 'sex' in pca_tab.columns:
                
                    good_male = [x for x in good_cells if x in pca_tab.loc[pca_tab.sex == 'male'].index]
                    good_female = [x for x in good_cells if x not in good_male]

                    plt.hist(subset[0].loc[event, good_male].dropna(), range=[0,1], density=True, 
                             color='skyblue', alpha=0.5, label='male')
                    plt.hist(subset[0].loc[event, good_female].dropna(), range=[0,1], density=True, 
                             color='forestgreen', alpha=0.5, label='female')
                    plt.title(event + ', ' + cell_type[j], fontsize=24)
                    plt.xlabel(r'$\Psi$', fontsize=24)
                    plt.ylabel('frequency', fontsize=24)
                    plt.legend(frameon=False)
                    plt.show()


                plt.figure()
                plt.scatter(pca_tab.loc[all_cells, 'PC1'], pca_tab.loc[all_cells, 'PC2'], c='grey', alpha=0.15)
                sc = plt.scatter(pca_tab.loc[good_cells, 'PC1'], pca_tab.loc[good_cells, 'PC2'], 
                        c=subset[0].loc[event, good_cells], vmin = 0, vmax=1, s=50)

                plt.colorbar(sc)

                plt.xlabel('PC1', fontsize=24)
                plt.ylabel('PC2', fontsize=24)
                plt.title(event + ', ' + cell_type[j], fontsize=24)

                plt.show()


                plt.scatter(np.log10(mrna_per_event.loc[event, good_cells]+1), 
                            subset[0].loc[event, good_cells], alpha=0.5)
                plt.xlabel('log10(mRNA counts)', fontsize=24)
                plt.ylabel(r'$\Psi$', fontsize=24)
                plt.title(event + ' ' + cell_type[j], fontsize=24)
                plt.show()


                plt.scatter(np.log10(read_counts.loc[event, good_cells]+1), 
                            subset[0].loc[event, good_cells], alpha=0.5)
                plt.xlabel('log10(SJ read counts)', fontsize=24)
                plt.ylabel(r'$\Psi$', fontsize=24)
                plt.title(event + ', ' + cell_type[j], fontsize=24)
                plt.show()
            i += 1

        j += 1



def hartigan_multiple_test(subset_list, cell_type, pca_tab, mrna_per_event, 
                     read_counts, quantile = 0.25, psi_lim = 0.25):
    pvals = []

    for subset in subset_list:

        for event in subset[0].index:

            ht = hartigan_test(subset[0].loc[event].dropna())
            pvals.append(ht[1])


    bh_correction = multipletests(pvals, method='fdr_bh')

    i = 0
    j = 0

    for subset in subset_list:

        for event in subset[0].index:

            if bh_correction[0][i]:


                plt.hist(subset[0].loc[event].dropna(), range=[0,1], density=True)
                plt.title(event + ', ' + cell_type[j], fontsize=24)
                plt.xlabel(r'$\Psi$', fontsize=24)
                plt.ylabel('frequency', fontsize=24)
                plt.show()

                good_cells = subset[0].loc[event].dropna().index
                
                if 'sex' in pca_tab.columns:
                
                    good_male = [x for x in good_cells if x in pca_tab.loc[pca_tab.sex == 'male'].index]
                    good_female = [x for x in good_cells if x not in good_male]

                    plt.hist(subset[0].loc[event, good_male].dropna(), range=[0,1], density=True, color='skyblue', alpha=0.5)
                    plt.hist(subset[0].loc[event, good_female].dropna(), range=[0,1], density=True, color='forestgreen', alpha=0.5)
                    plt.title(event + ', ' + cell_type[j], fontsize=24)
                    plt.xlabel(r'$\Psi$', fontsize=24)
                    plt.ylabel('frequency', fontsize=24)
                    plt.show()


                sc = plt.scatter(pca_tab.loc[good_cells, 'PC1'], pca_tab.loc[good_cells, 'PC2'], 
                        c=subset[0].loc[event, good_cells], vmin = 0, vmax=1)

                plt.colorbar(sc)

                plt.xlabel('PC1', fontsize=24)
                plt.ylabel('PC2', fontsize=24)
                plt.title(event + ', ' + cell_type[j], fontsize=24)

                plt.show()


                plt.scatter(np.log10(mrna_per_event.loc[event, good_cells]+1), 
                            subset[0].loc[event, good_cells], alpha=0.5)
                plt.xlabel('log10(mRNA counts)', fontsize=24)
                plt.ylabel(r'$\Psi$', fontsize=24)
                plt.title(event + ' ' + cell_type[j], fontsize=24)
                plt.show()


                plt.scatter(np.log10(read_counts.loc[event, good_cells]+1), 
                            subset[0].loc[event, good_cells], alpha=0.5)
                plt.xlabel('log10(SJ read counts)', fontsize=24)
                plt.ylabel(r'$\Psi$', fontsize=24)
                plt.title(event + ', ' + cell_type[j], fontsize=24)
                plt.show()
            i += 1

        j += 1

        
    

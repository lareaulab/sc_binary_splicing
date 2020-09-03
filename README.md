# sc_binary_splicing
Dataset analysis, simulations and figures from https://www.biorxiv.org/content/10.1101/2019.12.19.883256v1

**Coverage-dependent bias creates the appearance of binary splicing in single cells**

Single cell RNA sequencing provides powerful insight into the factors that determine each cell’s unique identity, including variation in transcription and RNA splicing among diverse cell types. Previous studies led to the surprising observation that alternative splicing outcomes among single cells are highly variable and follow a bimodal pattern: a given cell consistently produces either one or the other isoform for a particular splicing choice, with few cells producing both isoforms. Here we show that this pattern arises almost entirely from technical limitations. We analyzed single cell alternative splicing in human and mouse single cell RNA-seq datasets, and modeled them with a probabilistic simulator. Our simulations show that low gene expression and low capture efficiency distort the observed distribution of isoforms in single cells. This gives the appearance of a binary isoform distribution, even when the underlying reality is consistent with more than one isoform per cell. We show that accounting for the true amount of information recovered can produce biologically meaningful measurements of splicing in single cells. 

Single cell datasets in this analysis:

Dataset | Organism | Biological process | No. cells | Accession
---- | ---- | ---- | ---- | ----
Chen | mouse | mES motor neuron differentiation | 488 | GSE74155
Lescroart | mouse | cardiomyogenesis | 598 | GSE100471
Trapnell | human | skeletal myogenesis | 314 | GSE52529
Song | human | iPS motor neuron differentiation | 206 | GSE85908
Fletcher | mouse | olfactory neurogenesis | 849 | GSE95601

**Alignment**

Example specifying STAR options: [**data_processing/STAR.txt**](data_processing/STAR.txt)

We aligned the reads of each dataset using STAR 2.5.3 with two-pass mode and index overhang length adapted to the read length of each dataset. We used the hg38 genome annotation for the human RNA-seq datasets, and the mm10 annotation for the mouse datasets. 

**Gene expression**

Example specifying RSEM options: [**data_processing/RSEM.txt**](data_processing/RSEM.txt)

Gene expression levels in transcripts per million (TPM) were calculated by running RSEM on the BAM files produced by the STAR alignment. 

**Alternative splicing measurement**

[**data_processing/make_SJ_tables/**](data_processing/make_SJ_tables/)

[**data_processing/make_read_counts/**](data_processing/make_read_counts/)

We ran rMATS 3.2.5 on bulk human and mouse RNA-seq datasets from cell types matching the scRNA-seq datasets to find all annotated cassette exon alternative splicing events in each cell type. Then we used the SJ.out.tab files obtained from the scRNA-seq STAR alignment to obtain the splice junction reads compatible with the list of cassette exons found by rMATS.

**Alternative splicing analysis**

The output from the above steps is saved in [**data/**](data/) so that the downstream analysis can be reproduced.

The analyses of real scRNA-seq data are documented in an iPython notebook:

[**results/data_analysis/dataset_analysis_Figures_1_and_3.ipynb**](results/data_analysis/dataset_analysis_Figures_1_and_3.ipynb)

and the analysis of simulations and theoretical calculations are documented here:

[**results/simulations/Simulations.ipynb**](results/simulations/Simulations.ipynb)

[**results/simulations/Simulations_nuisance_factors.ipynb**](results/simulations/Simulations_nuisance_factors.ipynb)

[**results/simulations/Theoretical_calculations.ipynb**](results/simulations/Theoretical_calculations.ipynb)

R scripts and results for the simulations are found in:

[**results/simulations/**](results/simulations/standard_bias)

The simulations in this paper use SymSim. A copy of SymSim used to run the simulations is found in:

[**results/simulations/SymSim/**](results/simulations/standard_bias/SymSim)

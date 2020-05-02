library(foreach)
library(parallel)
library(doParallel)
library('fAsianOptions')
require("knitr")
library("devtools")
#load_all("SymSim")

# Multiple simulations for each pcr cDNA/pcr fragment intersection
# One gene length set and gene expression/meta for each simulation
# cDNA PCR cycles of cDNA and fragments vary between simulations.

set.seed(1)

ngenes <- 500
ncells <- 300
psi <- 0.5

cores=detectCores()
cl <- makeCluster(30) # adjust for number of cores
registerDoParallel(cl)

#for (sim in 1:30) {
foreach(sim = 1:30, .packages = c('knitr', 'fAsianOptions', 'devtools')) %dopar% {

   load_all("SymSim")
    
    data(gene_len_pool)

    gene_len_pool_selected <- gene_len_pool[((gene_len_pool >= quantile(gene_len_pool, probs=0.1)) & (gene_len_pool <= quantile(gene_len_pool, probs=0.9)))]
    
    exon_lengths <- scan('/mnt/lareaulab/cfbuenabadn/RNASeq/Human/Botvinnik/splicing/exon_lengths.txt')
    
    load("SymSim/data/len2nfrag.RData")
    
    mrna <- SimulateTrueCounts(ncells_total=ncells, ngenes=ngenes, evf_type="one.population", 
                               Sigma=0.5, randseed=sim, prop_hge = 0.0)

    suffix1 <- paste('sim', sim, sep='_')

    counter <- 0
    len_excluded <- c()
    len_included <- c()
    while (counter < 500) {
        excl <- sample(gene_len_pool_selected, 1)
        exon_len <- sample(exon_lengths, 1)
        incl <- excl + exon_len

        if (as.character(incl) %in% rownames(len2nfrag)) {

            len_excluded <- c(len_excluded, excl)
            len_included <- c(len_included, incl)
            counter <- counter + 1   
        }
    }    

    gene_len <- c(len_included, len_excluded)
    
    lapply(gene_len, write, paste('sim5/gene_length', suffix1, 'txt', sep='.'), append=TRUE, ncolumns=1000)

    write.table(mrna$cell_meta, paste('sim5/meta_cell', suffix1, 'tab', sep='.'), 
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')

    I_matrix = floor((mrna$counts) / 2)
    E_matrix = floor((mrna$counts) / 2)

    mrna_matrix <- rbind(I_matrix, E_matrix)

  #rate_2PCR=0.8
  for (pcr1 in 10:25) {
  for (pcr2 in 5:15){

        suffix <- paste('sim', sim, 'pcr1', pcr1, 'pcr2', pcr2, sep='_')
      
        observed_reads <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=mrna$cell_meta, protocol="nonUMI", alpha_sd=0.002,                                                       gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.02, nPCR1=pcr1, nPCR2=pcr2)

        write.table(mrna_matrix, paste('sim5/true_counts', suffix, 'tab', sep='.'), 
                    quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')

        write.table(observed_reads[[1]], paste('sim5/observed_counts', suffix, 'tab', sep='.'), 
                    quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')

    print(pcr2)


  
  }   #print('pcr1')
      print(pcr1)
  
  }   #print('sim')
      print(sim)
}    
    

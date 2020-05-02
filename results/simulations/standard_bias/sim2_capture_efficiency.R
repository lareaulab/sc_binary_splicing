library(foreach)
library(parallel)
library(doParallel)
library('fAsianOptions')
require("knitr")
library("devtools")


# Multiple simulations for each capture efficiency value
# One gene length set and gene expression/meta for each simulation
# Capture efficiency varies in each simulation.

set.seed(1)

ngenes <- 500
ncells <- 300
psi <- 0.5

#sim_round <- 1

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
    
    
  for (capture in seq(from = 0.01, to = 0.25, by=0.01)){
    #foreach(capture = seq(from = 0.01, to = 0.25, by=0.01), .packages = c('knitr', 'fAsianOptions', 'devtools')) %dopar% {
      round_seed <- sim + (10000*capture) # just making sure it is reproducible, but rand seeds are all different
      
      suffix <- paste('sim', sim, 'cap', capture, sep='_')
      
      mrna <- SimulateTrueCounts(ncells_total=ncells, ngenes=ngenes, evf_type="one.population", 
                               Sigma=0.5, randseed=round_seed, prop_hge = 0.0)
      

#         gene_len_excluded <- sample(c(1:1000),500,replace=FALSE)+1000 #sample(gene_len_pool_200, 500, replace = TRUE)
#         gene_len_included <- gene_len_excluded + 150
#         gene_len <- c(gene_len_included, gene_len_excluded)
      
      
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

      

    
    lapply(gene_len, write, paste('sim2/gene_length', suffix, 'txt', sep='.'), append=TRUE, ncolumns=1000)

    write.table(mrna$cell_meta, paste('sim2/meta_cell', suffix, 'tab', sep='.'), 
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')

    I_matrix = floor((mrna$counts) / 2)
    E_matrix = floor((mrna$counts) / 2)

    mrna_matrix <- rbind(I_matrix, E_matrix)

    observed_reads <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=mrna$cell_meta, protocol="nonUMI", alpha_mean=capture, 
                                          alpha_sd=0.002, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.02)

    write.table(mrna_matrix, paste('sim2/true_counts', suffix, 'tab', sep='.'), 
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')

    write.table(observed_reads[[1]], paste('sim2/observed_counts', suffix, 'tab', sep='.'), 
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')

    print(capture)

  }
  print(sim)
}



    

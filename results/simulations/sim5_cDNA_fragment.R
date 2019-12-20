library('fAsianOptions')
require("knitr")
library("devtools")
load_all("SymSim")

# Multiple simulations for each pcr cDNA/pcr fragment intersection
# One gene length set and gene expression/meta for each simulation
# cDNA PCR cycles of cDNA and fragments vary between simulations.

ngenes <- 500
ncells <- 300
psi <- 0.5

for (sim in 1:50) {
    
    mrna <- SimulateTrueCounts(ncells_total=ncells, ngenes=ngenes, evf_type="one.population", 
                               Sigma=0.5, randseed=sim, prop_hge = 0.0)

    suffix1 <- paste('sim', sim, sep='_')

    gene_len <- sample(c(1:1000),1000,replace=FALSE)+1000
    
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
      
        observed_reads <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=mrna$cell_meta, protocol="nonUMI", alpha_sd=0.002,                                                       gene_len=gene_len, depth_mean=1e5, depth_sd=3e3, lenslope=0.01, nPCR1=pcr1, nPCR2=pcr2)

        write.table(mrna_matrix, paste('sim5/true_counts', suffix, 'tab', sep='.'), 
                    quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')

        write.table(observed_reads[[1]], paste('sim5/observed_counts', suffix, 'tab', sep='.'), 
                    quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')


  
  }   #print('pcr1')
      #print(pcr1)
  
  }   #print('sim')
      #print(sim)
}    
    

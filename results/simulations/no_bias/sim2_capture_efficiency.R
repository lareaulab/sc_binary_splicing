library('fAsianOptions')
require("knitr")
library("devtools")
load_all("SymSim")

# Multiple simulations for each capture efficiency value
# One gene length set and gene expression/meta for each simulation
# Capture efficiency varies in each simulation.

set.seed(1)

ngenes <- 500
ncells <- 300
psi <- 0.5

for (sim in 1:30) {
#for (sim in 26:50) {
    
    mrna <- SimulateTrueCounts(ncells_total=ncells, ngenes=ngenes, evf_type="one.population", 
                               Sigma=0.5, randseed=sim, prop_hge = 0.0)

    suffix1 <- paste('sim', sim, sep='_')

    gene_len <- sample(c(1:1000),1000,replace=FALSE)+1000
    
    lapply(gene_len, write, paste('sim2/gene_length', suffix1, 'txt', sep='.'), append=TRUE, ncolumns=1000)

    write.table(mrna$cell_meta, paste('sim2/meta_cell', suffix1, 'tab', sep='.'), 
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')

    I_matrix = floor((mrna$counts) / 2)
    E_matrix = floor((mrna$counts) / 2)

    mrna_matrix <- rbind(I_matrix, E_matrix)

  
  for (capture in seq(from = 0.01, to = 0.25, by=0.01)){

        suffix <- paste('sim', sim, 'cap', capture, sep='_')
      
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



    

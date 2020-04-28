library('fAsianOptions')
require("knitr")
library("devtools")
load_all("SymSim")

# Multiple simulations for each PSI value
# Using 1 to 500 molecules for each simulation
# One gene length set and meta for each simulation
# PSI varies in each simulation.

set.seed(1)
ngenes <- 500
ncells <- 300

mrna <- matrix(1, 500, 300)
vec <- 1:500
mrna <- sweep(mrna, MARGIN=1, vec, `*`)

for (sim in 1:30) {
    
  kon_evf <- sapply(1:10, function(x) paste('kon_evf', x, sep=''))
  koff_evf <- sapply(1:10, function(x) paste('koff_evf', x, sep=''))
  s_evf <- sapply(1:10, function(x) paste('s_evf', x, sep=''))
  meta_colnames <- c("cellid", "pop", kon_evf, koff_evf, s_evf)
                  
  suffix1 <- paste('sim', sim, sep='_')
  gene_len <- sample(c(1:1000),1000,replace=FALSE)+1000
  lapply(gene_len, write, paste('sim3/gene_length', suffix1, 'txt', sep='.'), append=TRUE, ncolumns=1000)
                  
  meta_cell <- matrix(1, 300, 32)
  meta_cell[,1] <- sapply(1:300, function(x) paste('cell',x, sep='_'))
  meta_cell <- data.frame(meta_cell)
  for (evf in 3:32) {meta_cell[,evf] <- rnorm(300, mean=1, sd=0.5)}
  colnames(meta_cell) <- meta_colnames
                  
  write.table(meta_cell, paste('sim3/meta_cell', suffix1, 'tab', sep='.'), 
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
  
  for (psi in seq(from = 0.05, to = 0.5, by=0.01)) {
    
    suffix <- paste('psi', psi, 'sim', sim, sep='_')
    
    I_matrix = matrix(nrow = 500, ncol = 300)
    E_matrix = matrix(nrow = 500, ncol = 300)
    psi_true_matrix = matrix(nrow = 500, ncol = 300)
    
    for (i in 1:500) {
      
      I <- mapply(function(x, y) rbinom(1, x, y), mrna[i,], psi)
      E <- mrna[i,] - I
      psi_hat <- I/mrna[i,]
      I_matrix[i,] = I
      E_matrix[i,] = E
      psi_true_matrix[i,] = psi_hat
      
    }
    
    mrna_matrix <- rbind(I_matrix, E_matrix)
    
    observed_reads <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=meta_cell, protocol="nonUMI", alpha_mean=0.1, 
                                          alpha_sd=0.002, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.02)
    
    
    write.table(mrna_matrix, paste('sim3/true_counts', suffix, 'tab', sep='.'), 
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
    
    write.table(observed_reads[[1]], paste('sim3/observed_counts', suffix, 'tab', sep='.'), 
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
    
    print(paste('psi', psi))
  }
  
}



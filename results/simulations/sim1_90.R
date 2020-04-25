library('fAsianOptions')
require("knitr")
library("devtools")
load_all("SymSim")

# Re-run
####

gene_len <- c(sample(c(1:1000),1000,replace=FALSE)+1000,
              sample(c(1:1000),1000,replace=FALSE)+1000,
              sample(c(1:1000),1000,replace=FALSE)+1000)

ngenes <- 1500
mrna <- SimulateTrueCounts(ncells_total=500, ngenes=ngenes, evf_type="one.population", Sigma=0.5, randseed=0, prop_hge = 0.0)

################################
# Alpha, Beta = unif, unimodal #
################################

I_matrix = matrix(nrow = 1500, ncol = 500)
E_matrix = matrix(nrow = 1500, ncol = 500)
psi_true_matrix = matrix(nrow = 1500, ncol = 500)

alpha_param = c(runif(500,5,30), runif(500,1,30), 1/runif(500,1,5))
beta_param = c(runif(500,5,30), 1/runif(500,1,5), runif(500,1,30))

for (i in 1:1500) {
  PSI <- rbeta(500, alpha_param[i], beta_param[i])
  I <- mapply(function(x, y) rbinom(1, x, y), mrna$counts[i,], PSI)
  E <- mrna$counts[i,] - I
  psi_hat <- I/mrna$counts[i,]
  I_matrix[i,] = I
  E_matrix[i,] = E
  psi_true_matrix[i,] = psi_hat
  
}

mrna_matrix <- rbind(I_matrix, E_matrix)

observed_reads_2 <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=mrna[[3]], protocol="nonUMI", alpha_mean=0.25,
                                        alpha_sd=0.00002, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.01)

observed_reads_1 <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=mrna[[3]], protocol="nonUMI", alpha_mean=0.1, 
                                        alpha_sd=0.00002, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.01)

observed_reads_5 <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=mrna[[3]], protocol="nonUMI", alpha_mean=0.5, 
                                         alpha_sd=0.00002, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.01)

observed_reads_9 <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=mrna[[3]], protocol="nonUMI", alpha_mean=0.9, 
                                         alpha_sd=0.00002, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.01)

observed_reads_99 <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=mrna[[3]], protocol="nonUMI", alpha_mean=0.99, 
                                         alpha_sd=0.00002, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.01)

observed_reads_100 <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=mrna[[3]], protocol="nonUMI", alpha_mean=1, 
                                        alpha_sd=0.00002, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.01)


write.table(mrna_matrix, '../../../sim1_extended/unimodal_true.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(observed_reads_2[[1]], '../../../sim1_extended/unimodal_observed_2.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(observed_reads_1[[1]], '../../../sim1_extended/unimodal_observed_1.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(observed_reads_5[[1]], '../../../sim1_extended/unimodal_observed_5.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(observed_reads_9[[1]], '../../../sim1_extended/unimodal_observed_9.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(observed_reads_99[[1]], '../../../sim1_extended/unimodal_observed_99.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(observed_reads_100[[1]], '../../../sim1_extended/unimodal_observed_100.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')


lapply(gene_len, write, "../../../sim1_extended/gene_length.txt", append=TRUE, ncolumns=1000)

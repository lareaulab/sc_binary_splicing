library('fAsianOptions')
require("knitr")
library("devtools")
load_all("SymSim")

####

set.seed(1)

data(gene_len_pool)

gene_len_pool_200 <- gene_len_pool[gene_len_pool > 200]
gene_len_pool_200 <- gene_len_pool_200[gene_len_pool_200 < 3500]

gene_len_excluded <- sample(gene_len_pool_200, 500, replace = FALSE)
gene_len_included <- gene_len_excluded + 150

gene_len <- c(gene_len_included, gene_len_excluded)

print(mean(gene_len))

ngenes <- 500
mrna <- SimulateTrueCounts(ncells_total=300, ngenes=ngenes, evf_type="one.population", Sigma=0.5, randseed=1, prop_hge = 0.0)
              
################################
# Alpha, Beta = unif, unimodal #
################################

I_matrix = floor((mrna$counts) / 2)
E_matrix = floor((mrna$counts) / 2)

mrna_matrix <- rbind(I_matrix, E_matrix)

print(dim(mrna_matrix))


observed_reads_9 <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=mrna[[3]], protocol="nonUMI", alpha_mean=0.90,
                                        alpha_sd=0.002, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.01,
                                        amp_bias_limit=c(-0.1, 0.1))
              
observed_reads_7 <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=mrna[[3]], protocol="nonUMI", alpha_mean=0.75,
                                        alpha_sd=0.002, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.01,
                                        amp_bias_limit=c(-0.1, 0.1))
              
observed_reads_5 <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=mrna[[3]], protocol="nonUMI", alpha_mean=0.5,
                                        alpha_sd=0.002, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.01,
                                        amp_bias_limit=c(-0.1, 0.1))
              
observed_reads_2 <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=mrna[[3]], protocol="nonUMI", alpha_mean=0.25,
                                        alpha_sd=0.002, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.01,
                                        amp_bias_limit=c(-0.1, 0.1))

observed_reads_1 <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=mrna[[3]], protocol="nonUMI", alpha_mean=0.1, 
                                        alpha_sd=0.002, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.01,
                                        amp_bias_limit=c(-0.1, 0.1))

observed_reads_05 <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=mrna[[3]], protocol="nonUMI", alpha_mean=0.05, 
                                         alpha_sd=0.002, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.01,
                                        amp_bias_limit=c(-0.1, 0.1))

observed_reads_02 <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=mrna[[3]], protocol="nonUMI", alpha_mean=0.02, 
                                         alpha_sd=0.002, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.01,
                                        amp_bias_limit=c(-0.1, 0.1))

observed_reads_01 <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=mrna[[3]], protocol="nonUMI", alpha_mean=0.01, 
                                         alpha_sd=0.002, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.01,
                                        amp_bias_limit=c(-0.1, 0.1))


write.table(mrna_matrix, 'sim0/unimodal_true.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
              
write.table(observed_reads_9[[1]], 'sim0/unimodal_observed_9.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(observed_reads_7[[1]], 'sim0/unimodal_observed_7.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(observed_reads_5[[1]], 'sim0/unimodal_observed_5.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(observed_reads_2[[1]], 'sim0/unimodal_observed_2.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(observed_reads_1[[1]], 'sim0/unimodal_observed_1.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(observed_reads_05[[1]], 'sim0/unimodal_observed_05.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(observed_reads_02[[1]], 'sim0/unimodal_observed_02.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(observed_reads_01[[1]], 'sim0/unimodal_observed_01.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
              
lapply(gene_len, write, "sim0/gene_length.txt", append=TRUE, ncolumns=1000)

library('fAsianOptions')
require("knitr")
library("devtools")
load_all("SymSim")

####

set.seed(1)


data(gene_len_pool)

gene_len_pool_selected <- gene_len_pool[((gene_len_pool >= quantile(gene_len_pool, probs=0.1)) & (gene_len_pool <= quantile(gene_len_pool, probs=0.9)))]

exon_lengths <- scan('/mnt/lareaulab/cfbuenabadn/RNASeq/Human/Botvinnik/splicing/exon_lengths.txt')

#gene_len_pool_selected <- gene_len_pool[gene_len_pool <= quantile(gene_len_pool, probs=0.9)]

#gene_len_pool_200 <- gene_len_pool[gene_len_pool > 200]
#gene_len_pool_200 <- gene_len_pool_200[gene_len_pool_200 < 3500]

# gene_len_excluded <- c(sample(gene_len_pool_200, 500, replace = FALSE),
#                        sample(gene_len_pool_200, 500, replace = FALSE),
#                        sample(gene_len_pool_200, 500, replace = FALSE))

# gene_len_included <- gene_len_excluded + 150


load("SymSim/data/len2nfrag.RData")

counter <- 0
len_excluded <- c()
len_included <- c()
while (counter < 1500) {
    excl <- sample(gene_len_pool_selected, 1)
    exon_len <- sample(exon_lengths, 1)
    incl <- excl + exon_len
    
    if (as.character(incl) %in% rownames(len2nfrag)) {
        
        len_excluded <- c(len_excluded, excl)
        len_included <- c(len_included, incl)
        counter <- counter + 1
        
    }
    
}    


# gene_len_excluded <- c(sample(gene_len_pool_selected, 500, replace = TRUE),
#                        sample(gene_len_pool_selected, 500, replace = TRUE),
#                        sample(gene_len_pool_selected, 500, replace = TRUE))

# exon_len_excluded <- c(sample(exon_lengths, 500, replace = TRUE),
#                        sample(exon_lengths, 500, replace = TRUE),
#                        sample(exon_lengths, 500, replace = TRUE))

# gene_len_included <- gene_len_excluded + exon_len_excluded

gene_len <- c(len_included, len_excluded)
print(mean(gene_len))

ngenes <- 1500
mrna <- SimulateTrueCounts(ncells_total=300, ngenes=ngenes, evf_type="one.population", Sigma=0.5, randseed=1, prop_hge = 0.0)
              
################################
# Alpha, Beta = unif, unimodal #
################################

I_matrix = matrix(nrow = 1500, ncol = 300)
E_matrix = matrix(nrow = 1500, ncol = 300)
psi_true_matrix = matrix(nrow = 1500, ncol = 300)
              
alpha_param = c(runif(500,5,30), runif(500,1,30), 1/runif(500,1,5))
beta_param = c(runif(500,5,30), 1/runif(500,1,5), runif(500,1,30))

for (i in 1:1500) {
  PSI <- rbeta(300, alpha_param[i], beta_param[i])
  I <- mapply(function(x, y) rbinom(1, x, y), mrna$counts[i,], PSI)
  E <- mrna$counts[i,] - I
  psi_hat <- I/mrna$counts[i,]
  I_matrix[i,] = I
  E_matrix[i,] = E
  psi_true_matrix[i,] = psi_hat
  
}

mrna_matrix <- rbind(I_matrix, E_matrix)

observed_reads_9 <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=mrna[[3]], protocol="nonUMI", alpha_mean=0.90,
                                        alpha_sd=0.002, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.02)
              
observed_reads_7 <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=mrna[[3]], protocol="nonUMI", alpha_mean=0.75,
                                        alpha_sd=0.002, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.02)
              
observed_reads_5 <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=mrna[[3]], protocol="nonUMI", alpha_mean=0.5,
                                        alpha_sd=0.002, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.02)
              
observed_reads_2 <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=mrna[[3]], protocol="nonUMI", alpha_mean=0.25,
                                        alpha_sd=0.002, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.02)

observed_reads_1 <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=mrna[[3]], protocol="nonUMI", alpha_mean=0.1, 
                                        alpha_sd=0.002, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.02)

observed_reads_05 <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=mrna[[3]], protocol="nonUMI", alpha_mean=0.05, 
                                         alpha_sd=0.002, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.02)

observed_reads_02 <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=mrna[[3]], protocol="nonUMI", alpha_mean=0.02, 
                                         alpha_sd=0.002, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.02)

observed_reads_01 <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=mrna[[3]], protocol="nonUMI", alpha_mean=0.01, 
                                         alpha_sd=0.002, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.02)


write.table(mrna_matrix, 'sim1/unimodal_true.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
              
write.table(observed_reads_9[[1]], 'sim1/unimodal_observed_9.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(observed_reads_7[[1]], 'sim1/unimodal_observed_7.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(observed_reads_5[[1]], 'sim1/unimodal_observed_5.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(observed_reads_2[[1]], 'sim1/unimodal_observed_2.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(observed_reads_1[[1]], 'sim1/unimodal_observed_1.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(observed_reads_05[[1]], 'sim1/unimodal_observed_05.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(observed_reads_02[[1]], 'sim1/unimodal_observed_02.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(observed_reads_01[[1]], 'sim1/unimodal_observed_01.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
              
              
###############################
# Alpha, Beta = unif, bimodal #
###############################

I_matrix = matrix(nrow = 1500, ncol = 300)
E_matrix = matrix(nrow = 1500, ncol = 300)
psi_true_matrix = matrix(nrow = 1500, ncol = 300)
              
alpha_param = c(1/runif(500,5,30), runif(500,1,30), 1/runif(500,1,5))
beta_param = c(1/runif(500,5,30), 1/runif(500,1,5), runif(500,1,30))

for (i in 1:1500) {
  PSI <- rbeta(300, alpha_param[i], beta_param[i])
  I <- mapply(function(x, y) rbinom(1, x, y), mrna$counts[i,], PSI)
  E <- mrna$counts[i,] - I
  psi_hat <- I/mrna$counts[i,]
  I_matrix[i,] = I
  E_matrix[i,] = E
  psi_true_matrix[i,] = psi_hat
  
}

mrna_matrix <- rbind(I_matrix, E_matrix)

observed_reads_2 <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=mrna[[3]], protocol="nonUMI", alpha_mean=0.25,
					alpha_sd=0.002, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.02)

observed_reads_1 <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=mrna[[3]], protocol="nonUMI", alpha_mean=0.1, 
                                        alpha_sd=0.002, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.02)

observed_reads_05 <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=mrna[[3]], protocol="nonUMI", alpha_mean=0.05, 
                                         alpha_sd=0.002, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.02)

observed_reads_02 <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=mrna[[3]], protocol="nonUMI", alpha_mean=0.02, 
                                         alpha_sd=0.002, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.02)

observed_reads_01 <- True2ObservedCounts(true_counts=mrna_matrix, meta_cell=mrna[[3]], protocol="nonUMI", alpha_mean=0.01, 
                                         alpha_sd=0.002, gene_len=gene_len, depth_mean=1e5, depth_sd=3e4, lenslope=0.02)


write.table(mrna_matrix, 'sim1/bimodal_true.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(observed_reads_2[[1]], 'sim1/bimodal_observed_2.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(observed_reads_1[[1]], 'sim1/bimodal_observed_1.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(observed_reads_05[[1]], 'sim1/bimodal_observed_05.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(observed_reads_02[[1]], 'sim1/bimodal_observed_02.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
write.table(observed_reads_01[[1]], 'sim1/bimodal_observed_01.tab', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')


###########################################################################

lapply(gene_len, write, "sim1/gene_length.txt", append=TRUE, ncolumns=1000)

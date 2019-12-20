library(ggplot2)
library(gplots)
library(scone)
library(biomaRt)
library(RColorBrewer)
library(pheatmap)

library(destiny, quietly = TRUE)

library(mclust, quietly = TRUE)

library(SummarizedExperiment)
library(SingleCellExperiment)
library(clusterExperiment)
library(slingshot)


library(ggplot2)
library(gplots)
library(scone)
library(biomaRt)
library(pheatmap)




cc <- c(brewer.pal(9, "Set1"))
par(mar=c(1,1,1,1))




# 1. Get the data and metadata.

counts <- read.table('../data/chen/chen.tpm.gene_symbols.tab', header = TRUE, row.names = 1)
meta <- read.table('../data/chen/meta.tab', sep = '\t',
                   header = TRUE, row.names = 1)

neurogenesis_types = c('ES2i', 'ES', 'Epi', 'Motor neuron')
neurogenesis_cells = row.names(meta[meta$cell_type %in% neurogenesis_types,])
neuro_male_cells = row.names(meta[(meta$cell_type %in% neurogenesis_types) & meta$sex=='male',])
neuro_female_cells = row.names(meta[(meta$cell_type %in% neurogenesis_types) & meta$sex=='female',])



single_cells <- rownames(meta)#[(meta$cells_in_well == 1) & (meta$library_protocol == 'Single-cell')]

sc_counts <- counts[,single_cells]
sc_meta <- meta[single_cells, ]

#sc_counts <- sc_counts[, sc_meta$MBases > quantile(sc_meta$MBases)['25%']]
#sc_meta <- sc_meta[sc_meta$MBases > quantile(sc_meta$MBases)['25%'], ]


#######################


filtered_genes <- rownames(sc_counts[(rowSums(sc_counts > 20) > 20),])
filtered <- sc_counts[filtered_genes,]
log_counts <- log1p(filtered)


color_cell = c()
for (i in sc_meta$cell_type){
  if (i=='ES2i'){color_cell<-c(color_cell, 'red')}
  else if(i=='ES'){color_cell<-c(color_cell, 'orange')}
  else if(i=='Epi'){color_cell<-c(color_cell, 'forestgreen')}
  else {color_cell<-c(color_cell, 'blue')}
}


color_batch = c()
for (i in sc_meta$batch){
  if (i=='2015-10-19'){color_batch<-c(color_batch, 'red')}
  else {color_batch<-c(color_batch, 'blue')}
}

color_mb = c()
for (i in sc_meta$MBases){
  if (i<='111'){color_mb<-c(color_mb, 'red')}
  else if(i<='168'){color_mb<-c(color_mb, 'orange')}
  else if(i<='241'){color_mb<-c(color_mb, 'forestgreen')}
  else {color_mb<-c(color_mb, 'blue')}
}



pca <- prcomp(t(log_counts))
plot(pca$x[,1], pca$x[,2], col = color_cell, type = 'p',
     pch=16, xlab='PC1', ylab='PC2')

plot(pca$x[,1], pca$x[,2], col = color_batch, type = 'p',
     pch=16, xlab='PC1', ylab='PC2')

plot(pca$x[,1], pca$x[,2], col = color_mb, type = 'p',
     pch=16, xlab='PC1', ylab='PC2')



qcNames = c("MBases", "rsem_counts")

qc <- as.matrix(meta[,qcNames])





library(scRNAseq)


# 2.4.1 HK dropouts





mouse_hk <- c('Actb', 'Actg1', 'Anxa2', 'B2m', 'Calm1', 'Calm2', 'Cox6b1', 'Cox7c', 'Cstb',
              'Ctsd', 'Dynll1', 'Eef1a1', 'Eef2', 'Eif4g2', 'Ezr', 'Fau', 'Fth1', 'Gnb2l1',
              'H3f3a', 'Hsp90ab1', 'Hspa8', 'Naca', 'Nfkbia', 'Nme2', 'Pabpc1', 'Prdx1', 'Ptma',
              'Rpl13', 'Rpl13a', 'Rpl29', 'Rpl32', 'Rpl35', 'Rpl36al', 'Rpl37', 'Rpl5', 'Rpl8',
              'Rplp1', 'Rps11', 'Rps12', 'Rps13', 'Rps14', 'Rps16', 'Rps24', 'Rps25', 'Rps27a',
              'Rps5', 'Rps9', 'Tuba1b', 'Ubc','Uqcrh','Ywhaz')

hk <- intersect(mouse_hk, rownames(filtered))


mouse_pos <- c('Nova1', 'Nova2', 'Otx2', 'Fgf5', 'Eomes', 'Cer1', 'Foxa2', 'Map2', 'Dcx',
               'Nefl', 'Nefm', 'Isl1', 'Prph')

mpc <- hk <- intersect(mouse_pos, rownames(filtered))


hk_obs = rowMeans(log10(as.matrix(filtered[hk,])+1))
drop_outs = as.matrix(filtered[hk,]) == 0


poscon_good_cells = intersect(rownames(filtered), mpc)
negcon_good_cells = intersect(rownames(filtered), hk)


qcNames = c("MBases", "rsem_counts")
goodQC <- as.matrix(meta[,qcNames])

##################

# 2.4.2 False negative curves
batch = factor(sc_meta$batch)





filtered <- as.matrix(filtered)

num_reads = quantile(filtered[filtered > 0])[4]
num_cells = ncol(filtered)
is_common = rowSums(filtered >= num_reads ) >= num_cells

bio = factor(sc_meta$cell_type)







poscon_good_cells = intersect(rownames(filtered), mpc)
negcon_good_cells = intersect(rownames(filtered), hk)





########################################





ppq = scale(goodQC[,apply(goodQC,2,sd) > 0],center = TRUE,scale = TRUE)



my_scone <- SconeExperiment(filtered,
                            qc=ppq, bio = bio,
                            negcon_ruv = rownames(filtered) %in% negcon_good_cells,
                            poscon = rownames(filtered) %in% poscon_good_cells
)

EFF_FN = function (ei)
{
  sums = colSums(ei > 0)
  eo = t(t(ei)*sums/mean(sums))
  return(eo)
}

## ----- Scaling Argument -----

scaling=list(none=identity, # Identity - do nothing
             
             eff = EFF_FN, # User-defined function
             
             sum = SUM_FN, # SCONE library wrappers...
             tmm = TMM_FN, 
             uq = UQ_FN,
             fq = FQT_FN,
             deseq = DESEQ_FN)

########################################

my_scone <- scone(my_scone,
                  scaling=scaling,
                  #imputation = imputation, impute_args = impute_args,
                  k_qc=2, k_ruv = 3,
                  adjust_bio="no",
                  run=FALSE)


apply(get_params(my_scone),2,unique)
is_screened = ((get_params(my_scone)$imputation_method == "expect") &
                 (get_params(my_scone)$scaling_method %in% c("none",
                                                             "eff")))
my_scone = select_methods(my_scone,
                          rownames(get_params(my_scone))[!is_screened ])



BiocParallel::register(
  BiocParallel::SerialParam()
) # Register BiocParallel Serial Execution

my_scone <- scone(my_scone,
                  scaling=scaling,
                  #imputation = imputation, impute_args = impute_args,
                  run=TRUE,
                  eval_kclust = 2:6,
                  return_norm = "in_memory",
                  zero = "postadjust")





# View Metric Scores
head(get_scores(my_scone))

# View Mean Score Rank
head(get_score_ranks(my_scone))

# Extract normalized data from top method
out_norm = get_normalized(my_scone,
                          method = rownames(get_params(my_scone))[1], log=TRUE)



pc_obj = prcomp(apply(t(get_scores(my_scone)),1,rank),
                center = TRUE,scale = FALSE)
bp_obj = biplot_color(pc_obj,y = -get_score_ranks(my_scone),expand = .6)



bp_obj = biplot_color(pc_obj,y = -get_score_ranks(my_scone),expand = .6)

points(t(bp_obj[1,]), pch = 1, col = "red", cex = 1)
points(t(bp_obj[1,]), pch = 1, col = "red", cex = 1.5)

points(t(bp_obj[42,]), pch = 1, col = "blue", cex = 1)
points(t(bp_obj[42,]), pch = 1, col = "blue", cex = 1.5)

#points(t(bp_obj[rownames(bp_obj) == "none,uq,qc_k=2,no_bio,no_batch",]),
#       pch = 1, col = "blue", cex = 1)
#points(t(bp_obj[rownames(bp_obj) == "none,uq,qc_k=2,no_bio,no_batch",]),
#       pch = 1, col = "blue", cex = 1.5)

arrows(bp_obj[42,][1],
       bp_obj[42,][2],
       bp_obj[1,][1],
       bp_obj[1,][2],
       lty = 2, lwd = 2)

good_cells_norm <- get_normalized(my_scone,1, log=TRUE)



pca_norm <- prcomp(t(as.matrix(good_cells_norm)))


plot(pca_norm$x[,1], pca_norm$x[,2], col = color_cell, type = 'p', 
     pch=16, xlab='PC1', ylab='PC2')






best_norm <- good_cells_norm

nrows <- dim(best_norm)[1]
ncols <- dim(best_norm)[2]
tpm <- filtered
normalized <- best_norm

fano_factor <- apply(best_norm,1,function(x){var(x)/mean(x)})
fano_order <- names(fano_factor[order(fano_factor, decreasing=TRUE)])
top_fano <- fano_order[1:500]


norm_top_fano <- best_norm[top_fano,]

# PCA

pca <- prcomp(t(as.matrix(norm_top_fano)), scale. = FALSE)
rd1 <- pca$x[,1:2]

plot(rd1, col = color_cell, pch=16, asp = 1)

#####



pca <- prcomp(t(as.matrix(best_norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]

plot(rd1, col = color_cell, pch=16, asp = 1)




neurogenesis <- SingleCellExperiment(assays=list(normalized=best_norm, tpm=filtered),colData=sc_meta)

#####


dm <- DiffusionMap(t(best_norm))
rd2 <- cbind(DC1 = dm$DC1, DC2 = dm$DC2)

plot(rd2, col = color_cell, pch=16, asp = 1)


reducedDims(neurogenesis) <- SimpleList(PCA = rd1, DiffMap = rd2)


cl1 <- Mclust(rd1)$classification
colData(neurogenesis)$GMM <- cl1

plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)




cl2 <- kmeans(rd1, centers = 4)$cluster
colData(neurogenesis)$kmeans <- cl2

colData(neurogenesis)$bio = bio

plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1)





sce <- slingshot(neurogenesis,  reducedDim = 'PCA',  start.clus = 'ES2i', end.clus = 'Motor neuron', 
                 clusterLabels = 'cell_type')





summary(sce$slingPseudotime_1)




colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plot(reducedDims(sce)$PCA, col = colors[cut(sce$slingPseudotime_1,breaks=100)], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2)







t <- sce$slingPseudotime_1
boxplot(t[bio=='ES2i'], t[bio=='ES'], t[bio=='Epi'], t[bio=='Motor neuron'])


pseudo_curve <- data.frame(attributes(SlingshotDataSet(sce_psi))$curves$curve1$s)
pseudo_curve$ord <- attributes(SlingshotDataSet(sce_psi))$curves$curve1$ord

write.table(pseudo_curve, '../data/chen/psi_pseudotime.tab', quote = FALSE, sep='\t')







#write.table(best_norm, 'chen.norm_counts.tab', sep='\t', quote=FALSE)
#lapply(top_fano, write, "chen.topFano.txt", append=TRUE, ncolumns=1000)

single_cells_filtered = colnames(best_norm)

meta_t <- meta[single_cells_filtered,]
meta_t$pseudotime <- t

#write.table(meta_t, 'chen.meta_pseudotime.tab', sep='\t', quote=FALSE)




#meta_old <- read.table('/mnt/c/Users/ferna/Desktop/SingleCell/data/chen/chen.meta_pseudotime.tab', sep = '\t',
#                   header = TRUE, row.names = 1)


pseudo_curve <- data.frame(attributes(SlingshotDataSet(sce))$curves$curve1$s)
#pseudo_curve$ord <- attributes(SlingshotDataSet(sce))$curves$curve1$ord
colnames(pseudo_curve) <- c('line_1', 'line_2')

lineage_df <- rd1
lineage_df <- cbind(lineage_df, pseudo_curve)
lineage_df$pseudotime <- t
lineage_df$cell_type <- bio

write.table(lineage_df, '../data/chen/chen.pca.tab', quote = FALSE, sep='\t')

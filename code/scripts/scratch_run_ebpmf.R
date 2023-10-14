library(dplyr)
library(tidyverse)
library(data.table)
library(bedr)
library(readr)
library(NNLM)
library(CountClust)
library(fastTopics)
library(FNN)
library(robustbase)
library(smashr)
library(ebpmf)
library(Matrix)

source('/project2/mstephens/cfbuenabadn/gtex-stm/code/scripts/sgom.R')

load_coverage_counts <- function(tissue, gene){
    file_name <- paste0("coverage/counts/", tissue, "/", gene, ".csv.gz")
    df <- read_csv(file_name) %>%
    column_to_rownames(var = "Sample_ID") 
    
    rownames(df) <- sub(paste0("\\.", tissue), "", rownames(df))
        
    return (df)
    
}

group_matrix <- function(df, k = 10){
  n_groups <- ncol(df) %/% k
  
  # Split dataframe into groups of k columns and sum them
  output <- lapply(seq_len(n_groups), function(i) {
    group_start_col <- k * (i - 1) + 1
    group_end_col <- min(k * i, ncol(df))
    
    group_sum <- rowSums(df[, group_start_col:group_end_col], na.rm = TRUE)
    
    # Create a new data frame with one column and set its name
    group_df <- data.frame(group_sum)
    colnames(group_df) <- colnames(df)[group_start_col]
    
    return(group_df)
  })
  
  # Combine the output into a single dataframe
  output_df <- do.call(cbind, output)
  
  return(output_df)
}


plot_structure <- function(EL, gene_name, kfactors, annotation_){
    indis = rownames(EL)
    tissue_label  <- c()
    
    colnames(EL) <- sprintf("factor%d",seq(1:kfactors))
    colores = RColorBrewer::brewer.pal(kfactors+1,  "Paired")
    
    if (kfactors==2){
        colores <- c('#1f77b4','#ff7f0e')
    }
    
   StructureGGplot(EL, annotation = annotation_,
                      palette = colores, figure_title = gene_name, yaxis_label='samples',
                      axis_tick = list(axis_ticks_length = 0.1, 
                                       axis_ticks_lwd_y = 1, 
                                       axis_ticks_lwd_x = 1, 
                                       axis_label_size = 12, axis_label_face = "bold"),
                         legend_title_size = 12, legend_key_size = 1, legend_text_size = 12)
}


args = commandArgs(trailingOnly=TRUE)
gene_name = args[1]
K = as.integer(args[2])
plot1 = args[3]
plot2 = args[4]
png1 = args[5]
png2 = args[6]
output = args[7]

tissues <- c('Brain_Anterior_cingulate_cortex_BA24', 
             'Brain_Cortex', 
             'Brain_Frontal_Cortex_BA9', 
             'Brain_Putamen_basal_ganglia', 
             'Heart_Atrial_Appendage', 
             'Liver', 
             'Lung', 
             'Muscle_Skeletal', 
             'Skin_Not_Sun_Exposed_Suprapubic', 
             'Whole_Blood')

df1 <- load_coverage_counts('Brain_Anterior_cingulate_cortex_BA24', gene_name)
df2 <- load_coverage_counts('Brain_Cortex', gene_name)
df3 <- load_coverage_counts('Brain_Frontal_Cortex_BA9', gene_name)
df4 <- load_coverage_counts('Brain_Putamen_basal_ganglia', gene_name)
df5 <- load_coverage_counts('Heart_Atrial_Appendage', gene_name)
df6 <- load_coverage_counts('Liver', gene_name)
df7 <- load_coverage_counts('Lung', gene_name)
df8 <- load_coverage_counts('Muscle_Skeletal', gene_name)
df9 <- load_coverage_counts('Skin_Not_Sun_Exposed_Suprapubic', gene_name)
df10 <- load_coverage_counts('Whole_Blood', gene_name)

counts <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10)



samples <- read_tsv('config/samples.tsv') %>%
    column_to_rownames(var = "X1") %>%
    filter((group=='train') & (tissue_id %in% tissues)) 

counts <- counts %>% filter(rownames(counts) %in% rownames(samples))
nsamples <- dim(counts)[1]
nbases <- dim(counts)[2]

TotalRowCounts <- counts %>% rowSums()
rowCoverage <- (TotalRowCounts/nbases)

counts <- counts[(TotalRowCounts > 1000),]


if (nbases > 40000){
  scaling <- (nbases/20000) %>% round()
  scaling <- min(scaling, 50)
  counts <- counts %>% group_matrix(scaling)
}

nsamples_filtered <- dim(counts)[1]


counts <- cbind(rpois(nsamples_filtered, 0.1), counts, rpois(nsamples_filtered, 0.1)) %>% as.matrix()

fit_ebpmf = ebpmf_identity(counts,K=K,tol = 1e-3,maxiter = 100, init = 'fasttopics')

fit_sgom_nugget = cluster.mix(counts,K=K,tol=1e-3,maxit = 100,nugget=TRUE)
init <- list(F_init=t(fit_sgom_nugget$phi), L_init=fit_sgom_nugget$pi)

fit_ebpmf_sgom_init = ebpmf_identity(counts,K=K,tol = 1e-3,maxiter = 100, init = init)

indis <- counts %>% rownames() 
annotation_ = data.frame(sample_id = indis,
                         tissue_label = factor(samples[indis,'tissue_id']))

pdf(file=plot1)
plot_structure(fit_ebpmf$EL, gene_name, K, annotation_)
dev.off()

png(file=png1)
plot_structure(fit_ebpmf$EL, gene_name, K, annotation_)
dev.off()

pdf(file=plot2)
plot_structure(fit_ebpmf_sgom_init$EL, gene_name, K, annotation_)
dev.off()

png(file=png2)
plot_structure(fit_ebpmf_sgom_init$EL, gene_name, K, annotation_)
dev.off()


saveRDS(list(region=gene_name,
             annotation = annotation_,
             fit_ebpmf_EF = fit_ebpmf$EF,
             fit_ebpmf_EF_smooth = fit_ebpmf$EF_smooth,
             fit_ebpmf_EL = fit_ebpmf$EL,
             fit_sgom_nugget_phi = fit_sgom_nugget$phi,
             fit_sgom_nugget_pi = fit_sgom_nugget$pi,
             fit_ebpmf_sgom_init_EF = fit_ebpmf_sgom_init$EF,
             fit_ebpmf_sgom_init_EF_smooth = fit_ebpmf_sgom_init$EF_smooth,
             fit_ebpmf_sgom_init_EL = fit_ebpmf_sgom_init$EL,
             coords = colnames(counts),
             samples = rownames(counts)
            ),
        file=output
       )
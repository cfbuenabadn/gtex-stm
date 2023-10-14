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
# K = as.integer(args[2])
# plot1 = args[3]
# plot2 = args[4]
png1 = args[2]
png2 = args[3]
png3 = args[4]
png4 = args[5]
# png5 = args[6]
# png6 = args[7]
# png7 = args[8]
# png8 = args[9]

test_png1 = args[6]#[10]
test_png2 = args[7]#[11]
test_png3 = args[8]#[12]
test_png4 = args[9]#[13]
# test_png5 = args[14]
# test_png6 = args[15]
# test_png7 = args[16]
# test_png8 = args[17]

output = args[10]#[18]

print(output)

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

counts <- counts[(TotalRowCounts > 100),]


if (nbases > 35000){
  scaling <- (nbases/20000) %>% round()
  scaling <- min(scaling, 100)
  counts <- counts %>% group_matrix(scaling)
}

nsamples_filtered <- dim(counts)[1]

##########################################################################
# Selecting test vs train samples

select_test_tissue <- function(samples, counts, tissue, K=5) {
    female_samples <- samples %>% filter(rownames(samples) %in% rownames(counts)) %>%
       filter(tissue_id == tissue, sex == 'female') %>% rownames()
    female_test_samples <- female_samples[1:K]
    
    male_samples <- samples %>% filter(rownames(samples) %in% rownames(counts)) %>%
       filter(tissue_id == tissue, sex == 'male') %>% rownames()
    male_test_samples <- male_samples[1:K]
    
    test_samples <- c(female_test_samples, male_test_samples)
    return(test_samples)
}

get_test_samples <- function(samples, counts, tissues, K=5){
    test_samples <- c()
    for (tissue in tissues){
        tissue_test_samples <- select_test_tissue(samples, counts, tissue)
        test_samples <- c(test_samples, tissue_test_samples)
    }
    return (test_samples)
}

test_samples <- get_test_samples(samples, counts, tissues)
train_samples <- rownames(counts)[!(rownames(counts) %in% test_samples)]


##########################################################################

logCounts <- (counts+1) %>% log2() %>% round()

coords <- colnames(counts)
matrix_cols <- c('pois1', coords, 'pois2')

counts <- cbind(rpois(nsamples_filtered, 0.1), counts, rpois(nsamples_filtered, 0.1)) %>% as.matrix()
logCounts <- cbind(rpois(nsamples_filtered, 0.1), logCounts, rpois(nsamples_filtered, 0.1)) %>% as.matrix()

colnames(counts) <- matrix_cols
colnames(logCounts) <- matrix_cols


train_model <- function(counts, K, init='fasttopics'){
    
    if (init == 'sgom'){
        fit_sgom_nugget = cluster.mix(counts,K=K,tol=1e-3,maxit = 100,nugget=TRUE)
        init <- list(F_init=t(fit_sgom_nugget$phi), L_init=fit_sgom_nugget$pi)

    }
    fit_ebpmf = ebpmf_identity(counts,K=K,tol = 1e-3,maxiter = 100, init = init)
    return(fit_ebpmf)
}

fit_ebpmf_K2 <- train_model(counts[train_samples,], 2)
fit_ebpmf_sgom_K2 <- train_model(counts[train_samples,], 2, 'sgom')

fit_ebpmf_K3 <- train_model(counts[train_samples,], 3)
fit_ebpmf_sgom_K3 <- train_model(counts[train_samples,], 3, 'sgom')


#############################
    
# fit_ebpmf_K5 <- train_model(counts[train_samples,], 5)
# fit_ebpmf_sgom_K5 <- train_model(counts[train_samples,], 5, 'sgom')

# log_fit_ebpmf_K5 <- train_model(logCounts[train_samples,], 5)
# log_fit_ebpmf_sgom_K5 <- train_model(logCounts[train_samples,], 5, 'sgom')

############################

    
# indis <- counts %>% rownames() 
annotation_ = data.frame(sample_id = train_samples,
                         tissue_label = factor(samples[train_samples,'tissue_id']))

# pdf(file=plot1)
# plot_structure(fit_ebpmf$EL, gene_name, K, annotation_)
# dev.off()

png(file=png1)
plot_structure(fit_ebpmf_K2$EL, gene_name, 2, annotation_)
dev.off()

png(file=png2)
plot_structure(fit_ebpmf_sgom_K2$EL, gene_name, 2, annotation_)
dev.off()

png(file=png3)
plot_structure(fit_ebpmf_K3$EL, gene_name, 3, annotation_)
dev.off()

png(file=png4)
plot_structure(fit_ebpmf_sgom_K3$EL, gene_name, 3, annotation_)
dev.off()

# png(file=png5)
# plot_structure(fit_ebpmf_K5$EL, gene_name, 5, annotation_)
# dev.off()

# png(file=png6)
# plot_structure(fit_ebpmf_sgom_K5$EL, gene_name, 5, annotation_)
# dev.off()

# png(file=png7)
# plot_structure(log_fit_ebpmf_K5$EL, gene_name, 5, annotation_)
# dev.off()

# png(file=png8)
# plot_structure(log_fit_ebpmf_sgom_K5$EL, gene_name, 5, annotation_)
# dev.off()

get_model_list <- function(ebpmf_model){
    out_list <- list(EF = ebpmf_model$EF,
                     EF_smooth = ebpmf_model$EF_smooth,
                     EL = ebpmf_model$EL)
    return(out_list)
}



predict_factors <- function(counts, EF){
    FF <- as.matrix(EF)
    FF %>% dim() %>% print()
    fit_init = init_poisson_nmf(counts, F = FF, init.method = 'random')
    out = fit_poisson_nmf(counts,fit0=fit_init,update.factors = NULL)
    return(out)
}



get_model_k_list <- function(ebpmf_model, ebpmf_sgom_model, counts, test_samples){
    
    ebpmf_model_list <- get_model_list(ebpmf_model)
    ebpmf_sgom_model_list <- get_model_list(ebpmf_sgom_model)
    
    ebpmf_test <- predict_factors(counts[test_samples,], ebpmf_model_list$EF)
    ebpmf_test_L = (ebpmf_test$L)/rowSums(ebpmf_test$L)
    
    ebpmf_sgom_test <- predict_factors(counts[test_samples,], ebpmf_sgom_model_list$EF)
    ebpmf_sgom_test_L = (ebpmf_sgom_test$L)/rowSums(ebpmf_sgom_test$L)
    
    ebpmf_K <- list(ebpmf_model = ebpmf_model_list,
                    ebpmf_sgom_model = ebpmf_sgom_model_list,
                    ebpmf_test_L = ebpmf_test_L,
                    ebpmf_sgom_test_L = ebpmf_sgom_test_L)
    return(ebpmf_K)
    }

ebpmf_K2 <- get_model_k_list(fit_ebpmf_K2, fit_ebpmf_sgom_K2, counts, test_samples)
ebpmf_K3 <- get_model_k_list(fit_ebpmf_K3, fit_ebpmf_sgom_K3, counts, test_samples)
# ebpmf_K5 <- get_model_k_list(fit_ebpmf_K5, fit_ebpmf_sgom_K5, counts, test_samples)
# log_ebpmf_K5 <- get_model_k_list(log_fit_ebpmf_K5, log_fit_ebpmf_sgom_K5, counts, test_samples)

annotation_test = data.frame(sample_id = test_samples,
                         tissue_label = factor(samples[test_samples,'tissue_id']))

png(file=test_png1)
plot_structure(ebpmf_K2$ebpmf_test_L, gene_name, 2, annotation_test)
dev.off()

png(file=test_png2)
plot_structure(ebpmf_K2$ebpmf_sgom_test_L, gene_name, 2, annotation_test)
dev.off()

png(file=test_png3)
plot_structure(ebpmf_K3$ebpmf_test_L, gene_name, 3, annotation_test)
dev.off()

png(file=test_png4)
plot_structure(ebpmf_K3$ebpmf_sgom_test_L, gene_name, 3, annotation_test)
dev.off()

# png(file=test_png5)
# plot_structure(ebpmf_K5$ebpmf_test_L, gene_name, 5, annotation_test)
# dev.off()

# png(file=test_png6)
# plot_structure(ebpmf_K5$ebpmf_sgom_test_L, gene_name, 5, annotation_test)
# dev.off()

# png(file=test_png7)
# plot_structure(log_ebpmf_K5$ebpmf_test_L, gene_name, 5, annotation_test)
# dev.off()

# png(file=test_png8)
# plot_structure(log_ebpmf_K5$ebpmf_sgom_test_L, gene_name, 5, annotation_test)
# dev.off()

saveRDS(list(region=gene_name,
             annotation = annotation_,
             ebpmf_K2 = ebpmf_K2,
             ebpmf_K3 = ebpmf_K3,
             #ebpmf_K5 = ebpmf_K5,
             #log_ebpmf_K5 = log_ebpmf_K5,
             coords = coords,
             train_samples = train_samples,
             test_samples = test_samples
            ),
        file=output
       )
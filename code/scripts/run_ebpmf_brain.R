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

args = commandArgs(trailingOnly=TRUE)
gene_name = args[1]
K = as.integer(args[2])
log_counts = as.logical(args[3])
out_rds = args[4]
out_EF = args[5]
out_EF_smooth = args[6]
out_EL = args[7]
out_EL_test = args[8]

sgom_out_EF = args[9]
sgom_out_EF_smooth = args[10]
sgom_out_EL = args[11]
sgom_out_EL_test = args[12]


tissues <- c('Brain_Frontal_Cortex_BA9', 
             'Brain_Putamen_basal_ganglia', 
             'Liver', 
             'Muscle_Skeletal', 
             'Skin_Not_Sun_Exposed_Suprapubic',
             'Whole_Blood')


counts <- paste0("coverage/count_tables/", gene_name, ".tab.gz") %>%
    read_tsv() %>%
    column_to_rownames(var = "Sample_ID") 

    
    



test_samples <- character(0)
train_samples <- character(0)

# Loop through the elements of my_vector
for (element in row.names(counts)) {
  # Use regex to extract the prefix (X) and suffix (Y)
  match <- regmatches(element, regexpr("\\.(test|train)$", element))
  if (length(match) > 0) {
    # Remove ".Y" and add to the appropriate vector
    if (match == ".test") {
      test_samples <- c(test_samples, sub("\\.test$", "", element))
    } else if (match == ".train") {
      train_samples <- c(train_samples, sub("\\.train$", "", element))
    }
  }
}

     
     


row.names(counts) <- counts %>% row.names() %>% sub("\\..*", "", .)

##########################################################################

if (log_counts){
    print('monitor if logging counts')
    counts <- (counts+1) %>% log2() %>% round()
}
    
coords <- colnames(counts)
matrix_cols <- c('pois1', coords, 'pois2')

nsamples_filtered <- dim(counts)[1]

counts <- cbind(rpois(nsamples_filtered, 0.1), counts, rpois(nsamples_filtered, 0.1)) %>% as.matrix()
# logCounts <- cbind(rpois(nsamples_filtered, 0.1), logCounts, rpois(nsamples_filtered, 0.1)) %>% as.matrix()

colnames(counts) <- matrix_cols
# colnames(logCounts) <- matrix_cols


train_model <- function(counts, K, init='fasttopics'){
    
    if (init == 'sgom'){
        fit_sgom_nugget = cluster.mix(counts,K=K,tol=1e-3,maxit = 100,nugget=TRUE)
        init <- list(F_init=t(fit_sgom_nugget$phi), L_init=fit_sgom_nugget$pi)

    }
    fit_ebpmf = ebpmf_identity(counts,K=K,tol = 1e-3,maxiter = 100, init = init)
    return(fit_ebpmf)
}
     

fit_ebpmf <- train_model(counts[train_samples,], K)
fit_ebpmf_sgom <- train_model(counts[train_samples,], K, 'sgom')

# log_fit_ebpmf <- train_model(logCounts[train_samples,], K)
# log_fit_ebpmf_sgom <- train_model(logCounts[train_samples,], K, 'sgom')

############################
 
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

ebpmf_model <- get_model_k_list(fit_ebpmf, fit_ebpmf_sgom, counts, test_samples)
# log_ebpmf_model <- get_model_k_list(log_fit_ebpmf, log_fit_ebpmf_sgom, logCounts, test_samples)
     
sample_names <- row.names(counts)

samples <- read_tsv('config/samples.tsv') %>%
    column_to_rownames(var = "X1") %>%
    filter((group=='train') & (tissue_id %in% tissues)) 

annotation = data.frame(sample_id = sample_names,
                         tissue_label = factor(samples[sample_names,'tissue_id']))


saveRDS(list(gene=gene_name,
             K = K,
             annotation = annotation,
             ebpmf_model = ebpmf_model,
             coords = coords,
             train_samples = train_samples,
             test_samples = test_samples
            ),
        file=out_rds
       )
     
EF <- ebpmf_model$ebpmf_model$EF %>% as.data.frame()
factor_names <- colnames(EF)
EF['coords'] <- EF %>% row.names()
EF <- EF[,c('coords', factor_names)]
EF %>% as.data.frame() %>% write_tsv(out_EF)
     
EF_smooth <- ebpmf_model$ebpmf_model$EF_smooth %>% as.data.frame()
factor_names <- colnames(EF_smooth)
EF_smooth['coords'] <- EF_smooth %>% row.names()
EF_smooth <- EF_smooth[,c('coords', factor_names)]
EF_smooth %>% as.data.frame() %>% write_tsv(out_EF_smooth)
     
EL <- ebpmf_model$ebpmf_model$EL %>% as.data.frame()
factor_names <- colnames(EL)
EL['Sample_ID'] <- EL %>% row.names()
EL <- EL[,c('Sample_ID', factor_names)]
EL %>% as.data.frame() %>% write_tsv(out_EL)
     
EL_test <- ebpmf_model$ebpmf_test_L %>% as.data.frame()
factor_names <- colnames(EL_test)
EL_test['Sample_ID'] <- EL_test %>% row.names()
EL_test <- EL_test[,c('Sample_ID', factor_names)]
EL_test %>% as.data.frame() %>% write_tsv(out_EL_test)
     
     

     
EF <- ebpmf_model$ebpmf_sgom_model$EF %>% as.data.frame()
factor_names <- colnames(EF)
EF['coords'] <- EF %>% row.names()
EF <- EF[,c('coords', factor_names)]
EF %>% as.data.frame() %>% write_tsv(sgom_out_EF)
     
EF_smooth <- ebpmf_model$ebpmf_sgom_model$EF_smooth %>% as.data.frame()
factor_names <- colnames(EF_smooth)
EF_smooth['coords'] <- EF_smooth %>% row.names()
EF_smooth <- EF_smooth[,c('coords', factor_names)]
EF_smooth %>% as.data.frame() %>% write_tsv(sgom_out_EF_smooth)
     
EL <- ebpmf_model$ebpmf_sgom_model$EL %>% as.data.frame()
factor_names <- colnames(EL)
EL['Sample_ID'] <- EL %>% row.names()
EL <- EL[,c('Sample_ID', factor_names)]
EL %>% as.data.frame() %>% write_tsv(sgom_out_EL)

     
EL_test <- ebpmf_model$ebpmf_sgom_test_L %>% as.data.frame()
factor_names <- colnames(EL_test)
EL_test['Sample_ID'] <- EL_test %>% row.names()
EL_test <- EL_test[,c('Sample_ID', factor_names)]
EL_test %>% as.data.frame() %>% write_tsv(sgom_out_EL_test)
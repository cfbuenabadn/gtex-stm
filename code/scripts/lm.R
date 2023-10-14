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
ebpmf_rds = args[1]
out_tab = args[2]

prueba <- readRDS(ebpmf_rds)

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

samples <- read_tsv('config/samples.tsv')  %>%
    column_to_rownames(var = "X1") %>%
    filter((group=='train') & (tissue_id %in% tissues)) 


get.samples <- function(samples, tissue, sample_set){
    subset_samples <- samples %>% 
      filter(tissue_id == tissue) %>% 
      filter(row.names(.) %in% sample_set) %>% 
      row.names()
    return(subset_samples)
    }


get.factors <- function(EL, samples, factor='factor1'){
    factors <- EL %>% 
      as.data.frame() %>% 
      filter(row.names(EL) %in% samples) %>% 
      pull(factor) %>% 
      as.numeric()
    return(factors)
}


lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}





pair.maker <- function(x) {
  pairs <- combn(x, 2)
    
  filtered_pairs <- apply(pairs, 2, function(pair) pair[1] != pair[2])
                          
  result <- matrix(pairs[, filtered_pairs], ncol = sum(filtered_pairs))
                          
  result <- t(result)
                          
  colnames(result) <- NULL
                          
  return(result)
}
                          

                          
                          

fit_lm <- function(samples, t1, t2, sample_list, EL){
    tissue_1 <- get.samples(samples, t1, sample_list)
    tissue_2 <- get.samples(samples, t2, sample_list)
    
    t1_y <- get.factors(EL, tissue_1)
    t2_y <- get.factors(EL, tissue_2)
    y <- c(t1_y, t2_y)

    t1_x <- numeric(length(tissue_1))
    t2_x <- rep(1, length(tissue_2))
    x <- c(t1_x, t2_x)
    

    fit_train <- lm(y~x)
    
    return(fit_train)
}
                          

tissue_pair <- combn(tissues, 2) %>% t()
                          
annotation <- prueba$annotation

                          
EL_train <- prueba$ebpmf_model$ebpmf_model$EL
Y <- prueba$train_samples
tissue_list <- annotation$tissue_label %>% unique()
EL_train <- EL_train %>% as.data.frame() %>% filter(row.names(EL_train) %in% Y)
colnames(EL_train) <- c('factor1', 'factor2')
                          
EL_test <- prueba$ebpmf_model$ebpmf_test_L
Y <- prueba$test_samples
tissue_list <- annotation$tissue_label %>% unique()
EL_test <- EL_test %>% as.data.frame() %>% filter(row.names(EL_test) %in% Y)
colnames(EL_test) <- c('factor1', 'factor2')

                          
sgom_EL_train <- prueba$ebpmf_model$ebpmf_sgom_model$EL
Y <- prueba$train_samples
tissue_list <- annotation$tissue_label %>% unique()
sgom_EL_train <- sgom_EL_train %>% as.data.frame() %>% filter(row.names(sgom_EL_train) %in% Y)
colnames(sgom_EL_train) <- c('factor1', 'factor2')

                          
sgom_EL_test <- prueba$ebpmf_model$ebpmf_sgom_test_L
Y <- prueba$test_samples
tissue_list <- annotation$tissue_label %>% unique()
sgom_EL_test <- sgom_EL_test %>% as.data.frame() %>% filter(row.names(sgom_EL_test) %in% Y)
colnames(sgom_EL_test) <- c('factor1', 'factor2')
                 
pair_list <- character(0)
pval_train_list <- c()
pval_test_list <- c()
sgom_pval_train_list <- c()
sgom_pval_test_list <- c()

                          

for (i in 1:dim(tissue_pair)[1]){
    t1 <- tissue_pair[i,][1]
    t2 <- tissue_pair[i,][2]

    fit_test <- fit_lm(samples, t1, t2, prueba$test_samples, EL_test)
    pval_test <- lmp(fit_test)
    
    fit_train <- fit_lm(samples, t1, t2, prueba$train_samples, EL_train)
    pval_train <- lmp(fit_train)
    
    
    sgom_fit_test <- fit_lm(samples, t1, t2, prueba$test_samples, sgom_EL_test)
    sgom_pval_test <- lmp(sgom_fit_test)
    
    sgom_fit_train <- fit_lm(samples, t1, t2, prueba$train_samples, sgom_EL_train)
    sgom_pval_train <- lmp(sgom_fit_train)
    
    pair <- paste(t1, t2, sep=':')
    pair_list <- c(pair_list, pair)
    pval_test_list <- c(pval_test_list, pval_test)
    pval_train_list <- c(pval_train_list, pval_train)
    
    sgom_pval_test_list <- c(sgom_pval_test_list, sgom_pval_test)
    sgom_pval_train_list <- c(sgom_pval_train_list, sgom_pval_train)
    
}
                          
                          
X = data.frame(pair_list, 
               pval_train = pval_train_list,
               pval_test = pval_test_list,
               sgom_pval_train = sgom_pval_train_list,
               sgom_pval_test = sgom_pval_test_list)
                          
X %>% as.data.frame() %>% write_tsv(out_tab)
                          
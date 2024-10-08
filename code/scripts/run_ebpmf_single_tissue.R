Sys.setenv(VROOM_CONNECTION_SIZE=131072 * 100)
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
library(ebpmf.alpha)
library(Matrix)

source('/project2/mstephens/cfbuenabadn/gtex-stm/code/scripts/process_factor_functions.R')

args = commandArgs(trailingOnly=TRUE)
gene_name = args[1]
strand = as.character(args[2])

print(strand)
print(strand == 'plus')
print(strand == 'minus')

#dsf(x)

set.seed(123)

predict_factors <- function(counts, EF){
    FF <- as.matrix(EF)
    FF %>% dim() %>% print()
    set.seed(123)
    fit_init = init_poisson_nmf(counts, F = FF, init.method = 'random')
    out = fit_poisson_nmf(counts,fit0=fit_init,update.factors = NULL)
    return(out)
}

run_ebpmf <- function(counts, K) {

    ebpmf_out <- tryCatch(
        {
        lib_size <- rowSums(counts) %>% as.integer()
        set.seed(123)
        fit = ebpmf.alpha::ebpmf_identity(as.matrix(counts),K, lib_size=lib_size)
    
        ebpmf_out <- list(EF = fit$EF,
                          EF_smooth = fit$EF_smooth,
                          EL = fit$EL,
                          elbo = fit$elbo,
                          coords = colnames(counts),
                          samples = rownames(counts)
                          )
        },
        error = function(e) {
            ebpmf_out <- NULL
            }
        )

    return(ebpmf_out)
    
    }


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

counts <- paste0("coverage/counts_filtered/", gene_name, ".csv.gz") %>%
    read_csv() %>%
    column_to_rownames(var = "Sample_ID") 



coords_to_sort <- counts %>% colnames()
chrom <- (coords_to_sort[1] %>% strsplit(., ':'))[[1]][1]

coords_to_sort <- sapply(strsplit(coords_to_sort, ":"), function(x) as.integer(x[2])) %>% sort()
sorted_coords <- paste0(chrom, ':', coords_to_sort)
counts <- counts[sorted_coords]
                         
################

row.names(counts) <- counts %>% row.names() %>% sub("\\..*", "", .)

coords <- colnames(counts)

counts <- counts[rowSums(counts) >= 100, ]

kept_samples <- (counts %>% dim())[1]

if (kept_samples < 50){
    out_rds <- paste0('ebpmf_models/single_tissue/RDS/', gene_name, '.rds')

    saveRDS(list(gene=gene_name,
                 coords = coords,
                 ebpmf_run = NULL
                ),
            file=out_rds
       )
} else {


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


# attributes <- read_tsv(paste0('coverage/counts_filtered_stats/', gene_name, '.stats'), col_names = c('trait', 'quant'))
#nbases_total <- attributes %>% filter(trait == 'total_length') %>% pull(quant) %>% as.numeric()

nbases_total <- dim(counts)[2] ### CHANGED

if ((nbases_total > 15000) & (nbases_total < 60000)){
    print('warning: matrix is too big; merging groups of 3 base pairs to reduce size.')
    counts <- counts %>% group_matrix(., 3)
} else if (nbases_total > 60000){
    print('warning: matrix is too big; merging groups of 9 base pairs to reduce size.')
    counts <- counts %>% group_matrix(., 9)
}

##########################################################################


samples <- read_tsv('config/samples.tsv', skip=1, col_names = c('X1', 'tissue_id', 'sex', 'group')) %>% # read_tsv('config/samples.tsv') %>%
    column_to_rownames(var = "X1") %>%
    filter((group=='train') & (tissue_id %in% tissues)) 


train_and_test_ebpmf <- function(counts, train_samples, K, coords, strand){

    coords <- colnames(counts)

    matrix_cols <- c('pois1', coords, 'pois2')
    nsamples <- dim(counts)[1]

    set.seed(123)
    counts <- cbind(rpois(nsamples, 0.2), counts, rpois(nsamples, 0.2)) %>% as.matrix()
    colnames(counts) <- matrix_cols
    
    
    counts_test <- counts %>% 
                     as.data.frame() %>%
                     filter(!(rownames(counts) %in% train_samples)) %>% 
                     as.matrix()
    
    counts <- counts %>% 
                as.data.frame() %>%
                filter(rownames(counts) %in% train_samples) %>% 
                as.matrix()

    if ((dim(counts)[1] > 30) & (dim(counts_test)[1] > 1)){
        train_fit <- run_ebpmf(counts, K)
        test_predict <- predict_factors(counts_test, train_fit$EF)
    } else {
        train_fit <- NULL
        test_predict <- NULL
    }
    

    # train_fit <- run_ebpmf(counts, K)
    # test_predict <- predict_factors(counts_test, train_fit$EF)

    train_samples_ <- rownames(counts)
    test_samples_ <- rownames(counts_test)
    
    out <- list(train_fit=train_fit, 
                test_predict=test_predict,
                train_samples = train_samples_, 
                test_samples = test_samples_, 
                coords = coords)
    return(out)
    }



samples_wb <- samples %>%
    filter(tissue_id == 'Whole_Blood')

samples_ba9 <- samples %>%
    filter((group=='train') & (tissue_id == 'Brain_Frontal_Cortex_BA9')) 

samples_ms <- samples %>%
    filter(tissue_id == 'Muscle_Skeletal')

samples_sk <- samples %>%
    filter(tissue_id == 'Skin_Not_Sun_Exposed_Suprapubic')

print('hola 2')
ebpmf_wb <- train_and_test_ebpmf(counts, rownames(samples_wb), 3, coords, strand)
print('hola 3')
ebpmf_ba9 <- train_and_test_ebpmf(counts, rownames(samples_ba9), 3, coords, strand)
print('hola 4')
ebpmf_ms <- train_and_test_ebpmf(counts, rownames(samples_ms), 3, coords, strand)
print('hola 5')
ebpmf_sk <- train_and_test_ebpmf(counts, rownames(samples_sk), 3, coords, strand)


nbases <- dim(counts)[2]

print('hola final')

out_rds <- paste0('ebpmf_models/single_tissue/RDS/', gene_name, '.rds')

train_samples = rownames(samples)
test_samples = rownames(counts)[!(rownames(counts) %in% train_samples)]
print('hola save')

saveRDS(list(gene=gene_name,
             ebpmf_wb = ebpmf_wb,
             ebpmf_ba9 = ebpmf_ba9,
             ebpmf_ms = ebpmf_ms,
             ebpmf_sk = ebpmf_sk,
             coords = coords,
             train_samples = train_samples,
             test_samples = test_samples,
             strand = strand),
        file=out_rds
       )
print('hola done')


}

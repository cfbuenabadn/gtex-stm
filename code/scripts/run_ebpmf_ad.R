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
library(stats)

Sys.setenv(VROOM_CONNECTION_SIZE=5000720)

args = commandArgs(trailingOnly=TRUE)
gene_name = args[1]
chunk =  as.character(args[2])

group_matrix <- function(df, k = 9){
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


metadata = read_tsv('gao_data/metadata/merged_metadata.tab.gz') %>% data.frame()
rownames(metadata) <- metadata$specimenID


counts_ad <- paste0('gao_data/counts/', gene_name, '.csv.gz') %>%
    read_csv() %>%
    column_to_rownames(var = "Sample_ID") 

counts <- counts_ad[(metadata %>% filter(batch == 'batch2', libraryPreparationMethod == 'Kapa', sequencingBatch == 'NYGC3') %>% 
                        pull(specimenID)), ]# %>% as.matrix()

coords <- colnames(counts)
samples <- rownames(counts)
nbases_total <- dim(counts)[2]

if ((nbases_total > 15000) & (nbases_total < 60000)){
    print('warning: matrix is too big; merging groups of 3 base pairs to reduce size.')
    counts <- counts %>% group_matrix(., 3)
} else if (nbases_total > 60000){
    print('warning: matrix is too big; merging groups of 9 base pairs to reduce size.')
    counts <- counts %>% group_matrix(., 9)
}

filtered_coords <- colnames(counts)


counts <- counts %>% as.matrix()
lib_size_ <- rowSums(counts) %>% as.integer()
counts <- counts[lib_size_ >= 100,]
print(dim(counts))
kept_samples <- rownames(counts)

run_snmf <- function(counts, K, metadata) {

    coords <- colnames(counts)

    matrix_cols <- c('pois1', coords, 'pois2')
    nsamples <- dim(counts)[1]

    set.seed(123)
    counts <- cbind(rpois(nsamples, 0.2), counts, rpois(nsamples, 0.2)) %>% as.matrix()
    colnames(counts) <- matrix_cols

    lib_size <- rowSums(counts) %>% as.integer()
    #ad_samples <- rownames(counts)
    
    fit <- ebpmf_identity(counts, K, lib_size=lib_size)

    ad_samples <- rownames(counts)
    tissue_ad <- paste(metadata[ad_samples,'sequencingBatch'], 
                   metadata[ad_samples,'AD'], sep='.')
    
    tissues <- c(tissue_ad)#, tissue_gtex)
    EL <- fit$EL

    colnames(EL) <- paste0('factor', 1:K)

    annotation_ad = data.frame(sample_id = ad_samples,tissue_label = factor(tissue_ad))
    # annotation_gtex = data.frame(sample_id = gtex_samples,tissue_label = factor(tissue_gtex))
    
    # rownames(metadata) <- metadata$specimenID
    EL_annot <- data.frame(EL)
    EL_annot['AD_tissue_sex'] = annotation_ad$tissue_label
    EL_annot['AD'] = metadata[rownames(EL), 'AD']
    EL_annot['condition'] = metadata[rownames(EL), 'condition']
    EL_annot['tissue'] = metadata[rownames(EL), 'tissue_id']
    EL_annot['sex'] = metadata[rownames(EL), 'msex']
    EL_annot['counts'] = log1p(rowSums(counts[rownames(EL),]))

    out_list = list(ebpmf = fit, coords=coords)

    for (tissue in c('AC', 'DLPC', 'PCC')) {


        kw_factor1 <- kruskal.test(factor1 ~ condition, data = (EL_annot %>% filter(condition %in% c('no condition', 'AD'), tissue == tissue)))
        kw_factor2 <- kruskal.test(factor2 ~ condition, data = (EL_annot %>% filter(condition %in% c('no condition', 'AD'), tissue == tissue)))
        kw_factor3 <- kruskal.test(factor3 ~ condition, data = (EL_annot %>% filter(condition %in% c('no condition', 'AD'), tissue == tissue)))

        kw_factor1_sex0 <- kruskal.test(factor1 ~ condition, data = (EL_annot %>% filter(condition %in% c('no condition', 'AD'), tissue == tissue, sex==0)))
        kw_factor2_sex0 <- kruskal.test(factor2 ~ condition, data = (EL_annot %>% filter(condition %in% c('no condition', 'AD'), tissue == tissue, sex==0)))
        kw_factor3_sex0 <- kruskal.test(factor3 ~ condition, data = (EL_annot %>% filter(condition %in% c('no condition', 'AD'), tissue == tissue, sex==0)))

        kw_factor1_sex1 <- kruskal.test(factor1 ~ condition, data = (EL_annot %>% filter(condition %in% c('no condition', 'AD'), tissue == tissue, sex==1)))
        kw_factor2_sex1 <- kruskal.test(factor2 ~ condition, data = (EL_annot %>% filter(condition %in% c('no condition', 'AD'), tissue == tissue, sex==1)))
        kw_factor3_sex1 <- kruskal.test(factor3 ~ condition, data = (EL_annot %>% filter(condition %in% c('no condition', 'AD'), tissue == tissue, sex==1)))

        kruskal_tests <- list(
            kw_factor1 = kw_factor1,
            kw_factor2 = kw_factor2,
            kw_factor3 = kw_factor3
        )

        kruskal_tests_sex0 <- list(
            kw_factor1 = kw_factor1_sex0,
            kw_factor2 = kw_factor2_sex0,
            kw_factor3 = kw_factor3_sex0
        )

        kruskal_tests_sex1 <- list(
                kw_factor1 = kw_factor1_sex1,
                kw_factor2 = kw_factor2_sex1,
                kw_factor3 = kw_factor3_sex1
            )

        if (K == 5){
            kw_factor4 <- kruskal.test(factor4 ~ condition, data = (EL_annot %>% filter(condition %in% c('no condition', 'AD'), tissue == tissue)))
            kw_factor4_sex0 <- kruskal.test(factor4 ~ condition, data = (EL_annot %>% filter(condition %in% c('no condition', 'AD'), tissue == tissue, sex==0)))
            kw_factor4_sex1 <- kruskal.test(factor4 ~ condition, data = (EL_annot %>% filter(condition %in% c('no condition', 'AD'), tissue == tissue, sex==1)))

            kw_factor5 <- kruskal.test(factor5 ~ condition, data = (EL_annot %>% filter(condition %in% c('no condition', 'AD'), tissue == tissue)))
            kw_factor5_sex0 <- kruskal.test(factor5 ~ condition, data = (EL_annot %>% filter(condition %in% c('no condition', 'AD'), tissue == tissue, sex==0)))
            kw_factor5_sex1 <- kruskal.test(factor5 ~ condition, data = (EL_annot %>% filter(condition %in% c('no condition', 'AD'), tissue == tissue, sex==1)))

            kruskal_tests <- c(kruskal_tests,
                kw_factor4 = kw_factor4,
                kw_factor5 = kw_factor5
            )
            kruskal_tests_sex0 <- c(kruskal_tests_sex0,
                kw_factor4 = kw_factor4_sex0,
                kw_factor5 = kw_factor5_sex0
            )
            kruskal_tests_sex1 <- c(kruskal_tests_sex1,
                kw_factor4 = kw_factor4_sex1,
                kw_factor5 = kw_factor5_sex1
            )
        }

        tissue_kw <- list(kruskal_tests = kruskal_tests,
                         kruskal_tests_sex0 = kruskal_tests_sex0,
                         kruskal_tests_sex1 = kruskal_tests_sex1)

        out_list <- c(out_list, lst(!!tissue := tissue_kw))
        
    }

    return (out_list)

}

print(gene_name)


out_rds <- paste0('ebpmf_models/gao_models/RDS/',chunk, '/', gene_name, '.rds')

if (dim(counts) %>% is.null()){
    print('too few samples to run snmf')
    saveRDS(list(gene=gene_name,
             coords = coords,
             filtered_coords = filtered_coords,
             samples = samples,
             kept_samples = kept_samples
            ),
        file=out_rds
       )
} else if (dim(counts)[1] >= 30){
print('starting to run snmf')
snmf_3 <- run_snmf(counts, 3, metadata)
snmf_5 <- run_snmf(counts, 5, metadata)

print('finished running ebpmf')

saveRDS(list(gene=gene_name,
             snmf_3 = snmf_3,
             snmf_5 = snmf_5,
             coords = coords,
             filtered_coords = filtered_coords,
             samples = samples,
             kept_samples = kept_samples
            ),
        file=out_rds
       )

} else {
    print('too few samples to run snmf')
    saveRDS(list(gene=gene_name,
             coords = coords,
             filtered_coords = filtered_coords,
             samples = samples,
             kept_samples = kept_samples
            ),
        file=out_rds
       )
}

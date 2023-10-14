library(dplyr)
library(tidyverse)
library(data.table)


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

    
    
args = commandArgs(trailingOnly=TRUE)
gene_name = args[1]
out_name = args[2]
    
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

test_samples_names <- paste0(test_samples, ".test")
train_samples_names <- paste0(train_samples, ".train")
    
ordered_samples <- c(train_samples, test_samples)
sample_names <- c(train_samples_names, test_samples_names)
    
out <- rbind(counts[ordered_samples,])
# row.names(out) <- sample_names

coords <- colnames(out)
    
out['Sample_ID'] <- sample_names
    
out <- out[c('Sample_ID', coords)]
out %>% write_tsv(out_name)
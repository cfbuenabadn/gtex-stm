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

Sys.setenv(VROOM_CONNECTION_SIZE=50007200)

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


metadata = read_tsv('gregor_data/metadata/metadata_sberg.tab.gz') %>% data.frame()
rownames(metadata) <- metadata$participant_id_batch


counts <- paste0('gregor_data/sberger/counts/', gene_name, '.csv.gz') %>%
    read_csv() %>%
    column_to_rownames(var = "Sample_ID") 


# counts <- counts_ad[(metadata %>% filter(batch == 'batch2', libraryPreparationMethod == 'Kapa', sequencingBatch == 'NYGC3') %>% 
#                         pull(specimenID)), ]# %>% as.matrix()

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

seqmatic_counts <- counts[metadata %>% filter(within_site_batch_name == 'SeqMatic-1') %>% rownames(),]

invitae_samples <- c()
for (sample in metadata %>% filter(metadata$within_site_batch_name != 'SeqMatic-1') %>% rownames()){
    if (sample %in% rownames(counts)) {invitae_samples <- c(invitae_samples, sample)}
}

invitae_counts <- counts[invitae_samples,]


seqmatic_counts <- seqmatic_counts %>% as.matrix()
lib_size_ <- rowSums(seqmatic_counts) %>% as.integer()
seqmatic_counts <- seqmatic_counts[lib_size_ >= 100,]
print(dim(seqmatic_counts))
kept_samples_seqmatic <- rownames(seqmatic_counts)

invitae_counts <- invitae_counts %>% as.matrix()
lib_size_ <- rowSums(invitae_counts) %>% as.integer()
invitae_counts <- invitae_counts[lib_size_ >= 100,]
print(dim(invitae_counts))
kept_samples_invitae <- rownames(invitae_counts)

run_snmf <- function(counts, K) {

    coords <- colnames(counts)

    matrix_cols <- c('pois1', coords, 'pois2')
    nsamples <- dim(counts)[1]

    set.seed(123)
    counts <- cbind(rpois(nsamples, 0.2), counts, rpois(nsamples, 0.2)) %>% as.matrix()
    colnames(counts) <- matrix_cols

    lib_size <- rowSums(counts) %>% as.integer()
    #ad_samples <- rownames(counts)
    
    fit <- ebpmf_identity(counts, K, lib_size=lib_size)
    out_list = list(ebpmf = fit, coords=coords)

}

print(gene_name)


out_rds <- paste0('ebpmf_models/sberger_models/RDS/',chunk, '/', gene_name, '.rds')

if (dim(seqmatic_counts) %>% is.null()){
    print('too few samples to run snmf in seqmatic')
    seqmatic_object <- list(gene=gene_name,
             kept_samples = kept_samples_seqmatic
            )
} else if (dim(seqmatic_counts)[1] >= 30){
print('starting to run snmf in seqmatic')
snmf_3 <- run_snmf(seqmatic_counts, 3)
snmf_5 <- run_snmf(seqmatic_counts, 5)

print('finished running ebpmf in seqmatic')

seqmatic_object <- list(gene=gene_name,
             snmf_3 = snmf_3,
             snmf_5 = snmf_5,
             kept_samples = kept_samples_seqmatic
            )

} else {
    print('too few samples to run snmf in seqmatic')
    seqmatic_object <- list(gene=gene_name,
             kept_samples = kept_samples_seqmatic
            )
}


if (dim(invitae_counts) %>% is.null()){
    print('too few samples to run snmf in invitae')
    invitae_object <- list(gene=gene_name,
             kept_samples = kept_samples_invitae
            )
} else if (dim(invitae_counts)[1] >= 30){
print('starting to run snmf in invitae')
snmf_3 <- run_snmf(invitae_counts, 3)
snmf_5 <- run_snmf(invitae_counts, 5)

print('finished running ebpmf in invitae')

invitae_object <- list(gene=gene_name,
             snmf_3 = snmf_3,
             snmf_5 = snmf_5,
             kept_samples = kept_samples_invitae
            )

} else {
    print('too few samples to run snmf in invitae')
    invitae_object <- list(gene=gene_name,
             kept_samples = kept_samples_invitae
            )
}

print('wrapping up and saving')

saveRDS(list(gene=gene_name,
             coords = coords,
             filtered_coords = filtered_coords,
             samples = samples,
             invitae_run = invitae_object,
             seqmatic_run = seqmatic_object
            ),
        file=out_rds
       )

print('Done!')

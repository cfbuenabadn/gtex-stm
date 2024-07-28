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
library(ebpmf.alpha)#, lib = "/home/cnajar/R/x86_64-pc-linux-gnu-library/4.1-alt-packages/")
library(Matrix)
library(RNOmni)

Sys.setenv(VROOM_CONNECTION_SIZE=5000720)

source('/project2/mstephens/cfbuenabadn/gtex-stm/code/scripts/process_factor_functions.R')

gene_name = args[1]
out_rds = args[2]

metadata = read_tsv('gao_data/metadata/merged_metadata.tab.gz') %>% data.frame()
rownames(metadata) <- metadata$specimenID

set.seed(123)

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



counts_ad <- paste0('gao_data/counts/', gene_name, '.csv.gz') %>%
    read_csv() %>%
    column_to_rownames(var = "Sample_ID") 

coords <- colnames(counts_ad)
coords_start <- coords[1]
m <- length(coords)
coords_end <- coords[m] sub("^[^.]*\\:", "", .)

coords_range <- paste(coords_start, coords_end, sep='-')

bed_exons <- tabix(coords_range, 'Annotations/gencode.v44.primary_assembly.exons.bed.gz') %>%
    filter(gene_id == gene_name) #%>% filter(transcript_type %in% c('protein_coding'))

merged_bed <- bed_exons %>% bedr.merge.region()

intron_end <- (merged_bed %>% pull(start))[2:((merged_bed %>% dim())[1])] 
intron_start <- (merged_bed %>% pull(end))[1:((merged_bed %>% dim())[1]-1)]

skip_coords <- c()
for (i in 1:length(intron_start)){
    iend <- intron_end[i] 
    istart <- intron_start[i]
    intron_len <- iend - istart
    if (intron_len > 200) {
        skip_intron <- paste0('chr15:', (istart + 100):(iend - 100))
        skip_coords <- c(skip_coords, skip_intron)
    }
    }

keep_coords <- colnames(counts_ad)[!(colnames(counts_ad) %in% skip_coords)]
# counts_ad <- counts_ad[(metadata %>% filter(batch == 'batch1') %>% pull(specimenID)), keep_coords]# %>% as.matrix()

counts_ad <- counts_ad[(metadata %>% filter(batch == 'batch2', libraryPreparationMethod == 'Kapa') %>% pull(specimenID)), keep_coords]# %>% as.matrix()

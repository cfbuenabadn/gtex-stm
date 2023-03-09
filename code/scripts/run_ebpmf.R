library(dplyr)
library(tidyverse)
library(data.table)
library(ebpmf)
library(NNLM)

args = commandArgs(trailingOnly=TRUE)
tissue1 = args[1]
tissue2 = args[2]
geneName = args[3]
K = as.integer(args[4])
output = args[5]

dir <- '/project2/mstephens/cfbuenabadn/gtex-stm/code/'

samples <- read.table(paste0(dir, 'config/samples.tsv'))

train1_samples <- samples %>% 
                  filter(tissue_id == train1, group == 'train') %>% 
                  rownames() %>% 
                  paste0(train1, ".", .)

train2_samples <- samples %>% 
                  filter(tissue_id == train2, group == 'train') %>% 
                  rownames() %>% 
                  paste0(train2, ".", .)

tissue1_counts <- read.csv(paste0(dir, 'Counts/', tissue1, '/', geneName, '.Counts.csv.gz'),
                           row.names = 1, header=TRUE) %>% filter(rownames(.) %in% train1)
tissue2_counts <- read.csv(paste0(dir, 'Counts/', tissue2, '/', geneName, '.Counts.csv.gz'),
                           row.names = 1, header=TRUE) %>% filter(rownames(.) %in% train2)

X <- rbind(tissue1_counts, tissue2_counts) %>% as.matrix()

fit_ebpmf = ebpmf_identity(X,K=K)

saveRDS(list(gene=geneName,
             geneCounts=X,
             assays = c('RNASeq'),
             fit_ebpmf = fit_ebpmf
            ),
        file=paste(output,sep='')
       )





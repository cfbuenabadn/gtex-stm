#print(R.Version())
.libPaths( c( "/home/cnajar/R/x86_64-pc-linux-gnu-library/4.1" , .libPaths() ) )
print(.libPaths()) 

library(dplyr)
library(tidyverse)
library(data.table)
library(ebpmf, lib.loc="/project2/mstephens/cfbuenabadn/R.4.1")
library(NNLM, lib.loc="/project2/mstephens/cfbuenabadn/R.4.1")

args = commandArgs(trailingOnly=TRUE)
tissue_list = args[1]
geneName = args[2]
K = as.integer(args[3])
output = args[4]

tissue_list <- str_split(tissue_list, '[.]')[[1]]

dir <- '/project2/mstephens/cfbuenabadn/gtex-stm/code/'

samples <- read.table(paste0(dir, 'config/samples.tsv'))

X <- data.frame()

for (tissue in tissue_list){
    train_samples <- samples %>% 
                  filter(tissue_id == tissue, group == 'train') %>% 
                  rownames() %>% 
                  paste0(tissue, ".", .)
    tissue_counts <- read.csv(paste0(dir, 'Counts/', tissue, '/', geneName, '.Counts.csv.gz'),
                           row.names = 1, header=TRUE) %>% filter(rownames(.) %in% train_samples)
    X <- rbind(X, tissue_counts)
    }

X <- X %>% as.matrix()

fit_ebpmf = ebpmf_identity(X,K=K)

saveRDS(list(gene=geneName,
             geneCounts=X,
             assays = c('RNASeq'),
             tissues = tissue_list,
             fit_ebpmf = fit_ebpmf
            ),
        file=paste(output,sep='')
       )





#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : CalculatePi1_GetTraitPairs
# @created     : Thursday Jun 23, 2022 10:50:15 CDT
#
# @description :
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "10 ../code/scratch/PairwisePi1Traits.", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

Num_f_out_chunks <- as.numeric(args[1])
f_out_prefix <- args[2]
permutation_f_in <- args[-c(1:2)]

library(tidyverse)
# library(data.table)
library(qvalue)

PermutationPass.dat <- permutation_f_in %>%
  setNames(str_extract(., "(?<=QTLs/GTEx_10/)[^/]+")) %>%
  lapply(read_delim, delim=' ') %>% 
  bind_rows(.id="tissue_id") %>% 
  select(PC=tissue_id, phe_id, p_permutation=adj_beta_pval, beta=slope, beta_se=slope_se, singletrait_topvar=var_id, singletrait_topvar_chr = var_chr, singletrait_topvar_pos=var_from) %>% 
  #select(PC=tissue_id, grp_id, phe_id, p_permutation=adj_beta_pval, beta=slope, beta_se=slope_se, singletrait_topvar=var_id, singletrait_topvar_chr = var_chr, singletrait_topvar_pos=var_from) %>% 
  group_by(PC) %>% 
  mutate(FDR = qvalue(p_permutation)$qvalues) %>% 
  mutate(GeneLocus = str_extract(phe_id, "^[^.]+")) %>%
  #mutate(GeneLocus = grp_id) %>% 
  unite(Trait, PC, phe_id, sep=";")

print('Done reading tables')

dat.pairs <- PermutationPass.dat %>%
  filter(FDR<0.1) %>%
  left_join(., PermutationPass.dat, by="GeneLocus") %>%
  filter(!Trait.x==Trait.y) %>%
  separate(Trait.x, into=c("PC1", "P1"), sep=";") %>%
  separate(Trait.y, into=c("PC2", "P2"), sep=";") %>%
  group_by(GeneLocus) %>%
  mutate(Random = sample(Num_f_out_chunks,1))

print('Done creating data pairs')

dat.pairs %>%
  group_by(Random) %>%
  group_walk(~ write_tsv(.x, paste0(f_out_prefix, .y$Random, ".txt.gz")))

print('Done with everything')




library(bedr)
library(dplyr)
library(tidyverse)
library(qvalue)

args <- commandArgs(trailingOnly=TRUE)

PairwiseInput <- args[1]
Pi1Output <- args[2]

pi1_table <- read_tsv(PairwiseInput)
unique_combinations <- unique(pi1_table %>% select(PC1, PC2))

# # Initialize a vector to store the results
results <- vector("list", length = nrow(unique_combinations))

# # Iterate over each unique combination
for (i in 1:nrow(unique_combinations)) {
  pc1_value <- unique_combinations$PC1[i]
  pc2_value <- unique_combinations$PC2[i]

    print(c(pc1_value, pc2_value))
  
#   # Filter the data frame for the current combination
    tryCatch(
        {
            pair_pi1_1e1 <- 1-(pi1_table %>% filter(PC1 == pc1_value, PC2 == pc2_value, FDR.x <= 1e-1, trait.x.p.in.y < 1) %>%
     pull(trait.x.p.in.y) %>% pi0est())$pi0
            },
        error  = function(e){.GlobalEnv$pair_pi1_1e1 <- 1}
        )
  tryCatch(
        {
            pair_pi1_5e2 <- 1-(pi1_table %>% filter(PC1 == pc1_value, PC2 == pc2_value, FDR.x <= 5e-2, trait.x.p.in.y < 1) %>%
     pull(trait.x.p.in.y) %>% pi0est())$pi0
            },
        error  = function(e){.GlobalEnv$pair_pi1_5e2 <- 1}
        )

    tryCatch(
        {
            pair_pi1_1e2 <- 1-(pi1_table %>% filter(PC1 == pc1_value, PC2 == pc2_value, FDR.x <= 1e-2, trait.x.p.in.y < 1) %>%
     pull(trait.x.p.in.y) %>% pi0est())$pi0
            },
        error  = function(e){.GlobalEnv$pair_pi1_1e2 <- 1}
        )

    tryCatch(
        {
            print('hola')
            pair_pi1_1e4 <- 1-(pi1_table %>% filter(PC1 == pc1_value, PC2 == pc2_value, FDR.x <= 1e-4, trait.x.p.in.y < 1) %>%
     pull(trait.x.p.in.y) %>% pi0est())$pi0
            },
        error  = function(e){
            .GlobalEnv$pair_pi1_1e4 <- 1
        }
        )
    
  
  results[[i]] <- paste0(paste(pc1_value, pc2_value, pair_pi1_1e1, pair_pi1_5e2, pair_pi1_1e2, pair_pi1_1e4, sep = "\t"), '\n')
}

writeLines(unlist(results), Pi1Output)
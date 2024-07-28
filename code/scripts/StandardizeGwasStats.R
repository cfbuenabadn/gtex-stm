library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)

f_in <- args[1]
stats <- args[2]
f_out <- args[3]

chrom_list <- c(paste0('chr', 1:22), 'chrX')

if (stats == 'beta_se') {
    dat <- read_tsv(f_in, 
                col_names = c('chrom', 'start', 'end', 'var_id', 'var_alt_id', 'P', 'beta', 'SE')) %>%
                filter(chrom %in% chrom_list)
    beta <- dat$beta
    SE <- dat$SE   
} else if (stats == 'zscore') {
    dat <- read_tsv(f_in, 
                col_names = c('chrom', 'start', 'end', 'var_id', 'var_alt_id', 'P', 'Z', 'N', 'freq')) %>%
                filter(chrom %in% chrom_list)
    
    SE <- 1/sqrt((2*dat$freq)*(1-dat$freq)*(dat$N + ((dat$Z)^2)))
    beta <- dat$Z*SE
} else if (stats == 'beta_only') {
    dat <- read_tsv(f_in, 
                col_names = c('chrom', 'start', 'end', 'var_id', 'var_alt_id', 'P', 'beta', 'SE')) %>%
                filter(chrom %in% chrom_list)
    beta <- dat$beta
    Z <- (sign(beta))*abs(qnorm(dat$P/2))
    SE <- (dat$beta)/Z
} else if (stats == 'pvals_odd') {
    dat <- read_tsv(f_in, 
                col_names = c('chrom', 'start', 'end', 'var_id', 'var_alt_id', 'P', 'OR')) %>%
                filter(chrom %in% chrom_list)
    beta <- beta <- log(dat$OR)
    SE <- abs(beta/qnorm(dat$P/2))
} else if (stats == 'OR_SE') {
    dat <- read_tsv(f_in, 
                col_names = c('chrom', 'start', 'end', 'var_id', 'var_alt_id', 'P', 'OR', 'SE')) %>%
                filter(chrom %in% chrom_list)
    beta <- beta <- log(dat$OR)
    SE <- abs(beta/qnorm(dat$P/2))
} else {stop("Error in summary stats format")}



dat_out <- dat %>% 
    select(chrom, start, end, P, var_id, var_alt_id) %>%
    mutate(beta = beta) %>%
    mutate(SE=SE) %>%
    write_tsv(f_out)

library(sleuth)
library(dplyr)
library(readr)

args = commandArgs(trailingOnly=TRUE)
gp_file = args[1]
output = args[2]

s2c <- read_tsv(gp_file)
so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)
so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)

sleuth_table %>% write_tsv(output)
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
library(ebpmf)
library(Matrix)


plot_structure <- function(EL, gene_name, kfactors, annotation_){
    indis = rownames(EL)
    tissue_label  <- c()
    
    colnames(EL) <- sprintf("factor%d",seq(1:kfactors))
    colores = RColorBrewer::brewer.pal(kfactors+1,  "Paired")
    
    StructureGGplot(EL, annotation = annotation_,
                      palette = colores, figure_title = gene_name, yaxis_label='samples',
                      axis_tick = list(axis_ticks_length = 0.1, 
                                       axis_ticks_lwd_y = 1, 
                                       axis_ticks_lwd_x = 1, 
                                       axis_label_size = 12, axis_label_face = "bold"),
                         legend_title_size = 12, legend_key_size = 1, legend_text_size = 12)
}


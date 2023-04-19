#print(R.Version())
.libPaths( c( "/home/cnajar/R/x86_64-pc-linux-gnu-library/4.1" , .libPaths() ) )
print(.libPaths()) 

library(dplyr)
library(tidyverse)
library(data.table)
library(ebpmf)
library(NNLM)
library(Matrix)
library(CountClust)
library(fastTopics)
library(FNN)

args = commandArgs(trailingOnly=TRUE)
tissue_list = args[1]
geneName = args[2]
K = as.integer(args[3])
output = args[4]
structure_plot_tissue_id = args[5]
structure_plot_sex = args[6]
structure_plot_tissue_sex_id = args[7]

# Split the tissue_list string into a vector of tissue names
tissue_list <- str_split(tissue_list, '[.]')[[1]]

dir <- '/project2/mstephens/cfbuenabadn/gtex-stm/code/'
samples <- read.table(paste0(dir, 'config/samples.tsv'))

X <- data.frame()

# Loop over all tissues in tissue_list, and read the gene expression counts for each tissue
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

# Run the ebpmf algorithm on the gene expression data
set.seed(0)
fit_ebpmf = ebpmf_identity(X,K=K)

# Save model to RDS file
saveRDS(list(gene=geneName,
             geneCounts=X,
             assays = c('RNASeq'),
             tissues = tissue_list,
             fit_ebpmf = fit_ebpmf,
             coords = colnames(X),
             samples = rownames(X)
            ),
        file=paste(output,sep='')
       )


# Structure plot function
plot_structure <- function(fit, gene_name, kfactors, annotation = NULL, filter_by = 'tissue_id'){
    EL = fit$EL
    indis = rownames(EL)
    tissue_label  <- c()
    
    
    if (!is.null(annotation)) {
        sample_names <- EL %>% rownames() %>% sub("^[^.]+\\.", "", .)
        tissue_label <- annotation[sample_names, ] %>% pull(filter_by)  
    }
    
    annotation_ = data.frame(sample_id = indis,tissue_label = factor(tissue_label))
    print(head(annotation_))
    print(class(annotation_))
    print(class(annotation_$tissue_label)) 
    
    colnames(EL) <- sprintf("factor%d",seq(1:kfactors))
    colores = RColorBrewer::brewer.pal(kfactors,  "Paired")
    
   StructureGGplot(EL, annotation = annotation_,
                      palette = colores, figure_title = gene_name,
                      axis_tick = list(axis_ticks_length = 0.1, 
                                       axis_ticks_lwd_y = 1, 
                                       axis_ticks_lwd_x = 1, 
                                       axis_label_size = 12, axis_label_face = "bold"),
                         legend_title_size = 12, legend_key_size = 1, legend_text_size = 12)
}

options(repr.plot.width=7.5, repr.plot.height=10)


samples$tissue_sex_id <- paste(samples$tissue_id, samples$sex, sep = '-')

options(repr.plot.width=7.5, repr.plot.height=10)

# plot structure plot
pdf(file=structure_plot_tissue_id)
plot_structure(fit_ebpmf, geneName, K, annotation = samples, filter_by = 'tissue_id')
dev.off()

pdf(file=structure_plot_sex)
plot_structure(fit_ebpmf, geneName, K, annotation = samples, filter_by = 'sex')
dev.off()

pdf(file=structure_plot_tissue_sex_id)
plot_structure(fit_ebpmf, geneName, K, annotation = samples, filter_by = 'tissue_sex_id')
dev.off()


#print(R.Version())
.libPaths( c( "/home/cnajar/R/x86_64-pc-linux-gnu-library/4.1" , .libPaths() ) )
print(.libPaths()) 

library(dplyr)
library(tidyverse)
library(data.table)
library(ebpmf)
library(NNLM)

args = commandArgs(trailingOnly=TRUE)
tissue_list = args[1]
geneName = args[2]
K = as.integer(args[3])
output = args[4]
structure_plot = args[5]
factor_plot = args[6]

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
fit_ebpmf = ebpmf_identity(X,K=K)

# Save model to RDS file
saveRDS(list(gene=geneName,
             geneCounts=X,
             assays = c('RNASeq'),
             tissues = tissue_list,
             fit_ebpmf = fit_ebpmf
            ),
        file=paste(output,sep='')
       )


# Structure plot function
plot_structure <- function(fit, gene_name, kfactors){
    EL = fit$EL
    indis = rownames(EL)
    tissue_label  <- c()
    for (x in rownames(EL)){
        tissue_label <- c(tissue_label, strsplit(x, '[.]')[[1]][1])
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


# plot structure plot
png(filename=structure_plot)
plot_structure(fit_ebpmf, geneName, K)
dev.off()


# Plot factors

# Select height of the plot depending on the number of factors
getPlotHeight <- function(K) {
  if (K == 2) {
    output <- 6
  } else if (K == 3) {
    output <- 8
  } else if (K == 5) {
    output <- 11
  } else if (K == 10) {
    output <- 18
  } 
  return(output)
}

plot_height <- getPlotHeight(K)

# set color palette for factors
colores = RColorBrewer::brewer.pal(K,  "Paired")
options(repr.plot.width=15, repr.plot.height=plot_height)

png(filename=factor_plot)
par(mfrow=c(K, 1), mar=c(5, 5, 2, 1), oma=c(0, 0, 0, 0), cex.lab = 2, cex.axis = 1.5)
for(k in 1:K) {
  plot(runmed(fit_ebpmf$res$qf$Ef_smooth[,k],k=43),type='l',col=colores[k], lwd=2, 
     ylab= paste0("factor ", as.character(k)),
     xlab = 'position')
}
dev.off()




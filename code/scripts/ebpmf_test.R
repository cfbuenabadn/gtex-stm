.libPaths( c( "/home/cnajar/R/x86_64-pc-linux-gnu-library/4.1" , .libPaths() ) )
print(.libPaths())

library(dplyr)
library(tidyverse)
library(data.table)
library(ebpmf)#, lib.loc="/project2/mstephens/cfbuenabadn/R.4.1")
library(NNLM)#, lib.loc="/project2/mstephens/cfbuenabadn/R.4.1")
library(Matrix)
library(CountClust)
library(fastTopics)
library(FNN)

args = commandArgs(trailingOnly=TRUE)
input_rds = args[1]
tissue1 = args[2]
tissue2 = args[3]
geneName = args[4]
output = args[5]
# structure_plot = args[6]
# factor1_plot = args[7]
# factor2_plot = args[8]

K = 2

# plot_structure <- function(fit, gene_name, kfactors){
#     EL = fit$EL
#     indis = rownames(EL)
#     tissue_label  <- c()
#     for (x in rownames(EL)){
#         tissue_label <- c(tissue_label, strsplit(x, '[.]')[[1]][1])
#     }
#     annotation_ = data.frame(sample_id = indis,tissue_label = factor(tissue_label))
#     print(head(annotation_))
#     print(class(annotation_))
#     print(class(annotation_$tissue_label)) 
    
#     colnames(EL) <- sprintf("factor%d",seq(1:kfactors))
#     colores = RColorBrewer::brewer.pal(kfactors,  "Paired")
    
#    StructureGGplot(EL, annotation = annotation_,
#                       palette = colores, figure_title = gene_name,
#                       axis_tick = list(axis_ticks_length = 0.1, 
#                                        axis_ticks_lwd_y = 1, 
#                                        axis_ticks_lwd_x = 1, 
#                                        axis_label_size = 12, axis_label_face = "bold"),
#                          legend_title_size = 12, legend_key_size = 1, legend_text_size = 12)
# }

lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}

shrink_coords <- function(counts, coords) {
  counts_zero <- colSums(counts)
  K_left <- which.max(counts_zero > 0) %>% as.numeric()
  K_right <- (which.max(rev(counts_zero) > 0) - 1) %>% as.numeric()
  coords_for_EF <- coords[K_left:(length(coords) - K_right)]
  return(coords_for_EF)
}

dir <- '/project2/mstephens/cfbuenabadn/'
ebpmf_dir <- paste0(dir, 'gtex-stm/code/ebpmf_model/train_2tissues/')
annotation = read.table(paste0(dir, 'gtex-stm/code/config/samples.tsv'), sep='\t', header=1, row.names=1)

model = readRDS(input_rds)

coords <-  model$geneCounts %>% colnames()
new_coords <- shrink_coords(model$geneCounts, coords)

# png(filename=structure_plot)
# plot_structure(model$fit_ebpmf, geneName, K)
# dev.off()

# png(filename=factor1_plot)
# plot(model$fit_ebpmf$EF[,1],type='l',ylab='',main=paste(geneName, 'factor 1'))
# dev.off()

# png(filename=factor2_plot)
# plot(model$fit_ebpmf$EF[,2],type='l',ylab='',main=paste(geneName, 'factor 2'))
# dev.off()



dir <- '/project2/mstephens/cfbuenabadn/gtex-stm/code/'
samples <- read.table(paste0(dir, 'config/samples.tsv'))

tissue_list <- c(tissue1, tissue2)

EL <- model$fit_ebpmf$EL %>% as.data.frame()
tissue1_EL <- EL %>% as.data.frame() %>% filter(startsWith(rownames(.), tissue1))
tissue2_EL <- EL %>% as.data.frame() %>% filter(startsWith(rownames(.), tissue2))
KL_train <- KL.divergence(tissue1_EL[,'k1'], tissue2_EL[,'k1'], 1)

x <- EL %>% rownames() %>% startsWith(tissue1) %>% as.numeric()
y <- EL$k1 %>% as.numeric()

fit_train <- lm(y~x)
pval_train <- lmp(fit_train)

# female_test tissue1 vs tissue2

X <- data.frame()

for (tissue in tissue_list){
    train_samples <- samples %>%
                  filter(tissue_id == tissue, group == 'female_test') %>%
                  rownames() %>%
                  paste0(tissue, ".", .)
    tissue_counts <- read.csv(paste0(dir, 'Counts/', tissue, '/', geneName, '.Counts.csv.gz'),
                           row.names = 1, header=TRUE) %>% filter(rownames(.) %in% train_samples)
    X <- rbind(X, tissue_counts)
    }


X <- X %>% select(all_of(new_coords)) %>% as.matrix() 

FF <- as.matrix(model$fit_ebpmf$EF)
fit_init = init_poisson_nmf(X, F = FF, init.method = 'random')
out = fit_poisson_nmf(X,fit0=fit_init,update.factors = NULL)

EL_female_test <- data.frame(row.names = rownames(out$L))
EL_female_test['k1'] <- out$L[,'k1']/(out$L[,'k1'] + out$L[,'k2'])
EL_female_test['k2'] <- out$L[,'k2']/(out$L[,'k1'] + out$L[,'k2'])
rownames(EL_female_test) <- rownames(out$L)

tissue1_female_test_EL <- EL_female_test %>% as.data.frame() %>% filter(startsWith(rownames(.), tissue1))
tissue2_female_test_EL <- EL_female_test %>% as.data.frame() %>% filter(startsWith(rownames(.), tissue2))

KL_female_test <- KL.divergence(tissue1_female_test_EL[,'k1'], tissue2_female_test_EL[,'k1'], 1)

x <- EL_female_test %>% rownames() %>% startsWith(tissue1) %>% as.numeric()
y <- EL_female_test$k1 %>% as.numeric()
fit_female_test <- lm(y~x)
pval_female_test <- lmp(fit_female_test)

# female_test tissue1


female_test_list <- c('female_test', 'female_test2')

tissue_list <- c('Brain_Cortex', 'Muscle_Skeletal')
tissue <- 'Brain_Cortex'

X <- data.frame()

for (group_test in female_test_list){
    train_samples <- samples %>%
                  filter(tissue_id == tissue, group == group_test) %>%
                  rownames() %>%
                  paste0(tissue, ".", .)
    tissue_counts <- read.csv(paste0(dir, 'Counts/', tissue, '/', geneName, '.Counts.csv.gz'),
                           row.names = 1, header=TRUE) %>% filter(rownames(.) %in% train_samples)
    X <- rbind(X, tissue_counts)
    }

X <- X %>% select(all_of(new_coords)) %>% as.matrix() 

fit_init = init_poisson_nmf(X, F = FF, init.method = 'random')
out = fit_poisson_nmf(X,fit0=fit_init,update.factors = NULL)

EL_neg_tissue1 <- data.frame(row.names = rownames(out$L))
EL_neg_tissue1['k1'] <- out$L[,'k1']/(out$L[,'k1'] + out$L[,'k2'])
EL_neg_tissue1['k2'] <- out$L[,'k2']/(out$L[,'k1'] + out$L[,'k2'])
rownames(EL_neg_tissue1) <- rownames(out$L)

#x <- EL_neg %>% rownames() %>% . %in% female_test_samples %>% as.numeric()
y <- EL_neg_tissue1$k1 %>% as.numeric()
fit_neg_tissue1 <- lm(y~x)
pval_neg_tissue1 <- lmp(fit_neg_tissue1)



tissue <- 'Muscle_Skeletal'

X <- data.frame()

for (group_test in female_test_list){
    train_samples <- samples %>%
                  filter(tissue_id == tissue, group == group_test) %>%
                  rownames() %>%
                  paste0(tissue, ".", .)
    tissue_counts <- read.csv(paste0(dir, 'Counts/', tissue, '/', geneName, '.Counts.csv.gz'),
                           row.names = 1, header=TRUE) %>% filter(rownames(.) %in% train_samples)
    X <- rbind(X, tissue_counts)
    }

X <- X %>% select(all_of(new_coords)) %>% as.matrix() 

fit_init = init_poisson_nmf(X, F = FF, init.method = 'random')
out = fit_poisson_nmf(X,fit0=fit_init,update.factors = NULL)

EL_neg_tissue2 <- data.frame(row.names = rownames(out$L))
EL_neg_tissue2['k1'] <- out$L[,'k1']/(out$L[,'k1'] + out$L[,'k2'])
EL_neg_tissue2['k2'] <- out$L[,'k2']/(out$L[,'k1'] + out$L[,'k2'])
rownames(EL_neg_tissue2) <- rownames(out$L)

#x <- EL_neg %>% rownames() %>% . %in% female_test_samples %>% as.numeric()
y <- EL_neg_tissue2$k1 %>% as.numeric()
fit_neg_tissue2 <- lm(y~x)
pval_neg_tissue2 <- lmp(fit_neg_tissue2)

print('Is it when saving?')

saveRDS(list(gene=geneName,
             tissues = tissue_list,
             pval_train=pval_train,
             pval_female_test=pval_female_test,
             pval_neg_tissue1 = pval_neg_tissue1,
             pval_neg_tissue2 = pval_neg_tissue2
            ),
        file=paste(output, sep='')
       )


library(stm)
library(NNLM)

args = commandArgs(trailingOnly=TRUE)
geneCountsFile = args[1]
outputFile = args[2]

if (!dir.exists('stm_models')){
dir.create('stm_models')
}

geneName = strsplit(strsplit(geneCountsFile, "/")[[1]][2], ".Counts")[[1]][1]

geneCounts = read.csv(geneCountsFile, header=1, row.names=1)

fit_sgom = cluster.mix(geneCounts,K=5,tol=1e-3,maxit = 100,nugget=TRUE)

saveRDS(list(gene=geneName,
             geneCounts=geneCounts,
             assays = c('RNASeq'),
             fit_sgom = fit_sgom
            ),
        file=paste('stm_models/', geneName, '.sgom_K5.rds',sep='')
       )
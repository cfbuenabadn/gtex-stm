import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--gtf', type=str, required=True)
parser.add_argument('--clusters', type=str, required=True)
parser.add_argument('--output', type=str, required=True)

if __name__ == '__main__':
    
    args = parser.parse_args()
    gtf = args.gtf
    clusters = args.clusters
    output = args.output
    
    exons = pd.read_csv(gtf, sep='\t')

    clusters = pd.read_csv(clusters, sep=' ')
    
    print(clusters.head())
    
    clusters['chrom'] = clusters.index
    clusters = clusters.chrom

    cluster_bed = clusters.str.split(':', expand=True)
    cluster_bed.columns = ['chrom', 'start', 'end', 'clu']

    cluster_bed['strand'] = [x[-1] for x in cluster_bed.clu.str.split('_')]
    cluster_bed = cluster_bed[['chrom', 'start', 'end', 'strand', 'clu']]

    exons_str = exons.astype(str)
    cluster_str = cluster_bed.astype(str)

    X = pd.merge(cluster_str, exons_str,  how='left', 
             left_on=['chrom', 'end', 'strand'], 
             right_on = ['chr', 'start','strand']
            )

    Y = pd.merge(X, exons_str,  how='left', 
             left_on=['chrom', 'start_x', 'strand'], 
             right_on = ['chr', 'end','strand']
            )

    Z = Y.groupby('clu').gene_name_x.unique()
    Z2 = Y.groupby('clu').gene_name_y.unique()

    gene_overlap = [','.join(set([j for j in Z[i] if type(j)==str] + [j for j in Z2[i] if type(j)==str])) for i in range(len(Z))]

    idx_order = cluster_bed.clu.unique()
    chrom = cluster_bed.groupby('clu').chrom.first()
    start = cluster_bed.groupby('clu').start.min()
    end = cluster_bed.groupby('clu').end.max()
    strand = [x.split('_')[-1] for x in end.index]

    df = pd.DataFrame(index = end.index)
    df['chrom'] = chrom
    df['start'] = start
    df['end'] = end
    df['strand'] = strand
    df['gene'] = gene_overlap

    df = df.loc[idx_order]
    
    df.to_csv(output, sep='\t', index=True, header=True)
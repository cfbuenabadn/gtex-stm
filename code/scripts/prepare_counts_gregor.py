import pandas as pd
import gzip
import os
import numpy as np
import tabix
import sys
sys.path.append('/project2/mstephens/cfbuenabadn/gtex-stm/code/scripts')
from prepare_counts import *

def process_and_write_gregor(dataset, chrom, start, end, gene_name):
    output_name = f'gregor_data/{dataset}/counts/{gene_name}.csv.gz'
    bedgraph_dir = f'gregor_data/{dataset}/splicing_data/bedgraph/'
    file_count = 1
    with gzip.open(output_name, 'wt', compresslevel=6) as fh:
        
        file_list = sorted([bedgraph_dir + x for x in os.listdir(bedgraph_dir) if (x[-2:]=='gz')])
        
        for file_name in file_list:
            coord_str, counts_str = process_sample_gregor(dataset, file_name, chrom, start, end)
            if file_count == 1:
                fh.write(coord_str)
            fh.write(counts_str)
            file_count += 1

def process_sample_gregor(tissue, file_name, chrom, start, end):
    sample_name = file_name.split('/')[-1].split('.')[0]
    chrom_for_tabix = chrom
    gene_df = tabix_dataframe(file_name, chrom_for_tabix, start, end)
    coord_str, counts_str = get_counts_string_gregor(gene_df, sample_name, tissue, chrom, start, end)
    return coord_str, counts_str

def get_counts_string_gregor(gene_df, sample_name, tissue, chrom, start, end):
    coord_list = ['Sample_ID']
    sample_counts = [sample_name]
    for i, (index, row) in enumerate(gene_df.iterrows()):
        start_ = int(row.start)
        end_ = int(row.end)
        if i == 0:
            start_ = start
        if i == len(gene_df) - 1:
            end_ = end

        segment_length = end_ - start_
        segment_counts = int(float(row.counts))
        segment_to_add = [str(segment_counts)]*segment_length
        sample_counts += segment_to_add
        coord_list += list(range(start_, end_))

    coord_list = [chrom + ':' + str(x) if (x != 'Sample_ID') else x for x in coord_list]
    coord_str = ','.join(coord_list)+'\n'
    counts_str = ','.join(sample_counts)+'\n'
    
    return coord_str, counts_str
if __name__ == '__main__':

    arguments = sys.argv
    gene_name = arguments[1]
    dataset = arguments[2]
    
    selected_genes = pd.read_csv('../data/protein_coding_genes.bed.gz', sep='\t', #'../data/selected_genes.bed', sep='\t', 
                                 names = ['chrom', 'start', 'end', 'gene', 'gene_symbol', 'strand'])
    
    gene_bed = selected_genes.loc[selected_genes.gene==gene_name]
    print(gene_bed)
    chrom = str(gene_bed.chrom.iloc[0])
    start = int(gene_bed.start.iloc[0]) - 100
    end = int(gene_bed.end.iloc[0]) + 100

    process_and_write_gregor(dataset, chrom, start, end, gene_name)
    print('finished!')
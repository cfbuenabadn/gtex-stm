import pandas as pd
import gzip
import os
import numpy as np
import tabix
import sys
sys.path.append('/project2/mstephens/cfbuenabadn/gtex-stm/code/scripts')
from prepare_counts import *

def get_counts_string_gao_df(gene_df, chrom, start, end):
    coord_list = ['Sample_ID']
    sample_counts = []
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
    sample_counts = [int(x) for x in sample_counts]
    
    return coord_str, sample_counts

def get_counts_string_gao(file_prefix, chrom, start, end, strand):
    print(file_prefix)
    sample_name = file_prefix.split('.')[0]
    if strand == '+':
        file_1 = 'gao_data/bed_files/' + file_prefix + '.minus.bed.gz'
        file_2 = 'gao_data/bed_files/' + file_prefix + '.plus.second_read.bed.gz'
    elif strand == '-':
        file_1 = 'gao_data/bed_files/' + file_prefix + '.plus.bed.gz'
        file_2 = 'gao_data/bed_files/' + file_prefix + '.minus.second_read.bed.gz'
    else:
        raise Exception('strand error')

    print(file_1)
    print(file_2)

    gene_df_1 = tabix_dataframe(file_1, chrom, start, end)
    gene_df_2 = tabix_dataframe(file_2, chrom, start, end)

    coord_str_1, sample_counts_1 = get_counts_string_gao_df(gene_df_1, chrom, start, end)
    coord_str_2, sample_counts_2 = get_counts_string_gao_df(gene_df_2, chrom, start, end)

    assert coord_str_1 == coord_str_2

    sample_counts  = np.array(sample_counts_1) + np.array(sample_counts_2)
    sample_counts_str = [sample_name] + [str(x) for x in sample_counts]
    counts_str = ','.join(sample_counts_str)+'\n'

    return coord_str_1, counts_str

def process_and_write_ru(file_list, gene_series):

    gene_name = gene_series.gene
    chrom = gene_series.chrom
    start = int(gene_series.start) - 100
    end = int(gene_series.end) + 100
    strand  = gene_series.strand
    
    output_name = f'gao_data/counts/{gene_name}.csv.gz'
    bedgraph_dir = f'gao_data/bed_files/'
    file_count = 1

    with gzip.open(output_name, 'wb', compresslevel=6) as fh:
        
        #file_list = sorted([bedgraph_dir + x for x in os.listdir(bedgraph_dir) if (x[-2:]=='gz')])
        
        for file_name in file_list:
            coord_str, counts_str = get_counts_string_gao(file_name, chrom, start, end, strand)
            if file_count == 1:
                fh.write(coord_str.encode())
            fh.write(counts_str.encode())
            file_count += 1

if __name__ == '__main__':

    arguments = sys.argv
    gene_name = arguments[1]
    
    selected_genes = pd.read_csv('Annotations/gencode.v44.primary_assembly.protein_coding.genes.bed.gz', sep='\t', 
                                 names = ['chrom', 'start', 'end', 'gene', 'gene_symbol', 'strand'])

    gene_series = selected_genes.loc[selected_genes.gene==gene_name].iloc[0]
    
    first_read_plus = pd.Index([x.split('.plus.')[0] for x in os.listdir('gao_data/bed_files/') if x.endswith('.plus.bed.gz.tbi')])
    first_read_minus = pd.Index([x.split('.minus.')[0] for x in os.listdir('gao_data/bed_files/') if x.endswith('.minus.bed.gz.tbi')])
    
    second_read_plus = pd.Index([x.split('.plus.second_read.')[0] for x in os.listdir('gao_data/bed_files/') if x.endswith('.plus.second_read.bed.gz.tbi')])
    second_read_minus = pd.Index([x.split('.minus.second_read.')[0] for x in os.listdir('gao_data/bed_files/') if x.endswith('.minus.second_read.bed.gz.tbi')])

    file_intersection = first_read_plus.intersection(second_read_plus).intersection(first_read_minus).intersection(second_read_minus)

    print(len(file_intersection))

    process_and_write_ru(file_intersection, gene_series)
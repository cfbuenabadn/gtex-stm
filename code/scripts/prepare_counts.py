import pandas as pd
import gzip
import os
import numpy as np
import tabix
import sys

def process_and_write(tissue, chrom, start, end, gene_name):
    output_name = 'coverage/counts/' + tissue + '/' + gene_name + '.csv.gz'
    file_list = sorted(['coverage/bed/' + tissue + '/' + x for x in os.listdir('coverage/bed/' + tissue) if ('.tbi' not in x)])
    with gzip.open(output_name, 'wt', compresslevel=6) as fh:
        file_count = 1
        for file_name in file_list:
            coord_str, counts_str = process_sample(tissue, file_name, chrom, start, end)
            if file_count == 1:
                fh.write(coord_str)
            fh.write(counts_str)
            file_count += 1
    

def process_sample(tissue, file_name, chrom, start, end):
    sample_name = file_name.split('/')[-1].split('.')[0]
    chrom_for_tabix = chrom
    gene_df = tabix_dataframe(file_name, chrom_for_tabix, start, end)
    coord_str, counts_str = get_counts_string(gene_df, sample_name, tissue, chrom, start, end)
    return coord_str, counts_str


def tabix_dataframe(file_name, chrom, start, end):
    print(file_name)
    print(chrom, start, end)
    tb = tabix.open(file_name)
    print(type(chrom))
    print(type(start))
    print(type(end))
    print('just to make sure')
    print(file_name)
    print(chrom)
    print(start)
    print(end)
    records = tb.query(chrom, start, end)
    print('error is here')
    gene_df = pd.DataFrame(records, columns = ['chrom', 'start', 'end', 'counts'])
    return gene_df


def get_counts_string(gene_df, sample_name, tissue, chrom, start, end):
    coord_list = ['Sample_ID']
    sample_counts = [sample_name + '.' + tissue]
    for i, (index, row) in enumerate(gene_df.iterrows()):
        start_ = int(row.start)
        end_ = int(row.end)
        if i == 0:
            start_ = start
        if i == len(gene_df) - 1:
            end_ = end

        segment_length = end_ - start_
        segment_counts = int(row.counts)
        segment_to_add = [str(segment_counts)]*segment_length
        sample_counts += segment_to_add
        coord_list += list(range(start_, end_))

    coord_list = [chrom + ':' + str(x) if (x != 'Sample_ID') else x for x in coord_list]
    coord_str = ','.join(coord_list)+'\n'
    counts_str = ','.join(sample_counts)+'\n'
    
    return coord_str, counts_str

if __name__ == '__main__':
    #_, gene_name, tissue = 
    arguments = sys.argv
    gene_name = arguments[1]
    tissues = arguments[2:]
    
    selected_genes = pd.read_csv('../data/selected_genes.bed', sep='\t', 
                                 names = ['chrom', 'start', 'end', 'gene', 'gene_symbol', 'strand'])
    
    gene_bed = selected_genes.loc[selected_genes.gene==gene_name]
    print(gene_bed)
    chrom = str(gene_bed.chrom.iloc[0])
    start = int(gene_bed.start.iloc[0]) - 50
    end = int(gene_bed.end.iloc[0]) + 50

    for tissue in tissues:
        print('processing ' + tissue)
        process_and_write(tissue, chrom, start, end, gene_name)
        print('finished!')

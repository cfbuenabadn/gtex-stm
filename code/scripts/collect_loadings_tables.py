import numpy as np
import gzip
import os
from tqdm import tqdm
import rpy2
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
import pandas as pd

def load_ebpmf_gene(gene_id):
    readRDS = ro.r['readRDS']
    df = readRDS(f'ebpmf_models/filtered/RDS/{gene_id}.rds')
    with (ro.default_converter + pandas2ri.converter).context():
        pd_from_r_df = ro.conversion.get_conversion().rpy2py(df)
    return pd_from_r_df

def get_EL(rds, K):
    gene = str(rds['gene'][0])
    ebpmf_run = rds[f'ebpmf_{str(K)}']
    EL = pd.DataFrame(ebpmf_run['train_fit']['EL'])
    train_samples = list(ebpmf_run['train_samples']) 
    EL.index = train_samples
    EL.columns = [f'{gene}.factor_{i}' for i in range(1, K+1)]
    return EL

def process_gene(rds, samples_list, K, first_line = False):
    # rds = load_ebpmf_gene(gene)

    line_list = []
   
    if f'ebpmf_{str(K)}' in rds.keys():

        gene = str(rds['gene'][0])

        coords_list = list(rds['coords'])
        chrom, start = coords_list[0].split(':')
        
        _, end = coords_list[-1].split(':')
    
        strand_str = str(rds['strand'][0])
        if strand_str == 'plus':
            strand = '+'
        elif strand_str == 'minus':
            strand = '-'
        else:
            raise Exception('strand error')

        
        EL_df = get_EL(rds, K)
        EL_df.reindex(samples_list)

        if first_line:
            first_line = ['#Chr', 'start', 'end', 'pid', 'gid', 'strand'] + samples_list
            first_line = '\t'.join(first_line) + '\n'
            line_list.append(first_line.encode())

        for idx, row in EL_df.reindex(samples_list).T.iterrows():
            
            line = [chrom, start, end, idx, gene, strand] + [str(str(round(x, 4))) for x in row]
            line = ['NA' if (x == 'nan') else x for x in line]
            out_line = '\t'.join(line) + '\n'
            line_list.append(out_line.encode())
    return line_list

if __name__ == '__main__':
    genes = [x.split('.')[0] for x in os.listdir('ebpmf_models/filtered/RDS/')]

    rds = load_ebpmf_gene('ENSG00000112081')
    samples_list = list(rds['train_samples'])

    fh_10 = gzip.open('ebpmf_models/filtered/snmf_10/tables/EL.bed.gz', 'wb')
    fh_2 = gzip.open('ebpmf_models/filtered/snmf_2/tables/EL.bed.gz', 'wb')
    fh_3 = gzip.open('ebpmf_models/filtered/snmf_3/tables/EL.bed.gz', 'wb')
    fh_4 = gzip.open('ebpmf_models/filtered/snmf_4/tables/EL.bed.gz', 'wb')
    fh_5 = gzip.open('ebpmf_models/filtered/snmf_5/tables/EL.bed.gz', 'wb')
    
    first_line_2 = True
    first_line_3 = True
    first_line_4 = True
    first_line_5 = True
    first_line_10 = True
    
    for gene in tqdm(genes, position=0, leave=True):
        rds = load_ebpmf_gene(gene)
        row_list_2 = process_gene(rds, samples_list, 2, first_line_2)
        row_list_3 = process_gene(rds, samples_list, 3, first_line_3)
        row_list_4 = process_gene(rds, samples_list, 4, first_line_4)
        row_list_5 = process_gene(rds, samples_list, 5, first_line_5)
        row_list_10 = process_gene(rds, samples_list, 10, first_line_10)
        
        if len(row_list_2) > 0:
            fh_2.writelines(row_list_2)
            first_line_2 = False
            
        if len(row_list_3) > 0:
            fh_3.writelines(row_list_3)
            first_line_3 = False

        if len(row_list_4) > 0:
            fh_4.writelines(row_list_4)
            first_line_4 = False

        if len(row_list_5) > 0:
            fh_5.writelines(row_list_5)
            first_line_5 = False

        if len(row_list_10) > 0:
            fh_10.writelines(row_list_10)
            first_line_10 = False
            

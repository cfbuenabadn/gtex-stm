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
    df = readRDS(f'ebpmf_models/single_tissue/RDS/{gene_id}.rds')
    with (ro.default_converter + pandas2ri.converter).context():
        pd_from_r_df = ro.conversion.get_conversion().rpy2py(df)
    return pd_from_r_df

def get_EL(rds, tissue):
    gene = str(rds['gene'][0])
    ebpmf_run = rds[f'ebpmf_{tissue}']
    EL = pd.DataFrame(ebpmf_run['train_fit']['EL'])
    train_samples = list(ebpmf_run['train_samples']) 
    EL.index = train_samples
    K = EL.shape[1]
    EL.columns = [f'{gene}.factor_{i}' for i in range(1, K+1)]
    return EL

def process_gene(rds, samples_list, tissue, first_line = False):
    # rds = load_ebpmf_gene(gene)

    line_list = []
   
    if f'ebpmf_{tissue}' in rds.keys():

        try:

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
    
            
            EL_df = get_EL(rds, tissue)
            EL_df.reindex(samples_list)
    
            # print(EL_df.shape)
    
            if first_line:
                sample_list_line = ['-'.join(x.split('-')[:2]) for x in samples_list]
                first_line = ['#Chr', 'start', 'end', 'pid', 'gid', 'strand'] + sample_list_line
                first_line = '\t'.join(first_line) + '\n'
                line_list.append(first_line.encode())
    
            for idx, row in EL_df.reindex(samples_list).T.iterrows():
                
                line = [chrom, start, end, idx, gene, strand] + [str(str(round(x, 4))) for x in row]
                line = ['NA' if (x == 'nan') else x for x in line]
                # print(len(line))
                out_line = '\t'.join(line) + '\n'
                line_list.append(out_line.encode())
        except:
            line_list = []
    return line_list


if __name__ == '__main__':
    genes = [x.split('.')[0] for x in os.listdir('ebpmf_models/single_tissue/RDS/')]

    rds = load_ebpmf_gene('ENSG00000112081')
    samples_list_ba9 = list(rds['ebpmf_ba9']['train_samples'])
    samples_list_wb = list(rds['ebpmf_wb']['train_samples'])
    samples_list_ms = list(rds['ebpmf_ms']['train_samples'])
    samples_list_sk = list(rds['ebpmf_sk']['train_samples'])

    fh_ba9 = gzip.open('ebpmf_models/single_tissue/tables/BA9.EL.bed.gz', 'wb')
    fh_wb = gzip.open('ebpmf_models/single_tissue/tables/WB.EL.bed.gz', 'wb')
    fh_ms = gzip.open('ebpmf_models/single_tissue/tables/MS.EL.bed.gz', 'wb')
    fh_sk = gzip.open('ebpmf_models/single_tissue/tables/Skin.EL.bed.gz', 'wb')
    
    first_line_ba9 = True
    first_line_wb = True
    first_line_ms = True
    first_line_sk = True
    
    # with gzip.open('ebpmf_models/filtered/snmf_10/tables/EL.bed.gz', 'wb') as fh:
    for gene in tqdm(genes, leave=True, position=0):
        rds = load_ebpmf_gene(gene)
        row_list_ba9 = process_gene(rds, samples_list_ba9, 'ba9', first_line_ba9)
        row_list_wb = process_gene(rds, samples_list_wb, 'wb', first_line_wb)
        row_list_ms = process_gene(rds, samples_list_ms, 'ms', first_line_ms)
        row_list_sk = process_gene(rds, samples_list_sk, 'sk', first_line_sk)
        
        if len(row_list_ba9) > 0:
            fh_ba9.writelines(row_list_ba9)
            first_line_ba9 = False
            
        if len(row_list_wb) > 0:
            fh_wb.writelines(row_list_wb)
            first_line_wb = False

        if len(row_list_ms) > 0:
            fh_ms.writelines(row_list_ms)
            first_line_ms = False

        if len(row_list_sk) > 0:
            fh_sk.writelines(row_list_sk)
            first_line_sk = False
            

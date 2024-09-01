import numpy as np
import pandas as pd
import gzip
import os
import sys
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri

def load_ebpmf_gene(dataset, gene_id, n):
    readRDS = ro.r['readRDS']
    df = readRDS(f'/project2/mstephens/cfbuenabadn/gtex-stm/code/ebpmf_models/{dataset}/RDS/{str(n)}/{gene_id}.rds')
    with (ro.default_converter + pandas2ri.converter).context():
        pd_from_r_df = ro.conversion.get_conversion().rpy2py(df)

    output = pd_from_r_df

    return output
    
def get_EL(rds, k):
    col_names = [f'{gene}.factor_{str(i+1)}' for i in range(k)]
    EL = pd.DataFrame(rds[f'snmf_{str(k)}']['ebpmf']['EL'])
    EL.columns = col_names
    EL.index = list(rds['kept_samples'])

    return EL

def process_gene(dataset, gene_id, chunk, gene_annotation, k=3, first_line = False):
    line_list = []
    rds = load_ebpmf_gene(dataset, gene_id, chunk)
    samples = sorted(list(rds['samples']))

    if first_line:
        first_line = ['#Chr', 'start', 'end', 'pid', 'gid', 'strand'] + samples
        first_line = '\t'.join(first_line) + '\n'
        line_list.append(first_line.encode())

    gene_coords = gene_annotation.loc[gene_annotation.gene_id == gene]
    chrom = gene_coords['#Chr'].iloc[0]
    start = str(gene_coords['start'].iloc[0])
    end = str(gene_coords['end'].iloc[0])
    strand = str(gene_coords['strand'].iloc[0])

    
    
    if f'snmf_{str(k)}' in rds.keys():
        EL_3 = get_EL(rds, k)
        EL_3 = EL_3.reindex(samples).T


        for factor, row in EL_3.iterrows():
            line = [chrom, start, end, factor, gene_id + f'.{chunk}', strand] + [str(str(round(x, 4))) for x in row]
            line = ['NA' if (x == 'nan') else x for x in line]
            out_line = '\t'.join(line) + '\n'
            line_list.append(out_line.encode())
            
    return line_list
            

if __name__ == '__main__':

    # Process Gao data
    gao_chunks = os.listdir('/project2/mstephens/cfbuenabadn/gtex-stm/code/ebpmf_models/gao_models/RDS/')

    first_pass = True
    gene_annotation = pd.read_csv('Annotations/gencode.v44.primary_assembly.protein_coding.genes.bed.gz', sep='\t')

    with gzip.open('/project2/mstephens/cfbuenabadn/gtex-stm/code/ebpmf_models/gao_models/tables/snmf3_EL.tab.gz', 'wb') as fh_3:
        with gzip.open('/project2/mstephens/cfbuenabadn/gtex-stm/code/ebpmf_models/gao_models/tables/snmf5_EL.tab.gz', 'wb') as fh_5:

            for chunk in gao_chunks:
                genes = os.listdir('/project2/mstephens/cfbuenabadn/gtex-stm/code/ebpmf_models/gao_models/RDS/' + chunk)
                for gene_id in genes:
                    if gene_id.startswith('ENSG'):
                        gene = gene_id.split('.')[0]
                        snmf3_list = process_gene('gao_models', gene, chunk, gene_annotation, 3, first_pass)
                        snmf5_list = process_gene('gao_models', gene, chunk, gene_annotation, 5, first_pass)
                        fh_3.writelines(snmf3_list)
                        fh_5.writelines(snmf5_list)
                        first_pass = False


    # sberger_chunks = os.listdir('/project2/mstephens/cfbuenabadn/gtex-stm/code/ebpmf_models/sberger_models/RDS/')

    # first_pass = True
    # gene_annotation = pd.read_csv('Annotations/gencode.v44.primary_assembly.protein_coding.genes.bed.gz', sep='\t')

    # with gzip.open('/project2/mstephens/cfbuenabadn/gtex-stm/code/ebpmf_models/sberger_models/tables/snmf3_EL.tab.gz', 'wb') as fh_3:
    #     with gzip.open('/project2/mstephens/cfbuenabadn/gtex-stm/code/ebpmf_models/sberger_models/tables/snmf5_EL.tab.gz', 'wb') as fh_5:

    #         for chunk in sberger_chunks:
    #             genes = os.listdir('/project2/mstephens/cfbuenabadn/gtex-stm/code/ebpmf_models/sberger_models/RDS/' + chunk)
    #             for gene_id in genes:
    #                 if gene_id.startswith('ENSG'):
    #                     gene = gene_id.split('.')[0]
    #                     snmf3_list = process_gene('sberger_models', gene, chunk, gene_annotation, 3, first_pass)
    #                     snmf5_list = process_gene('sberger_models', gene, chunk, gene_annotation, 5, first_pass)
    #                     fh_3.writelines(snmf3_list)
    #                     fh_5.writelines(snmf5_list)
    #                     first_pass = False
    

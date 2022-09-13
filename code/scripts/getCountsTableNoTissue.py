import numpy as np
import pandas as pd
import argparse


def remove_tissue(data, annotation, tissue):
    idx_without_tissue = data.index[annotation.loc[data.index].tissue_id != tissue]
    data_without_tissue = data.loc[idx_without_tissue]
    return data_without_tissue


def one_sample_for_tissue(data, annotation, tissue, seed = 1):
    idx_tissue = data.index[annotation.loc[data.index].tissue_id == tissue]
    np.random.seed(seed)
    tissue_sample = np.random.choice(idx_tissue)
    idx_plus_one_tissue = list(data.index[annotation.loc[data.index].tissue_id != tissue])
    idx_plus_one_tissue.append(tissue_sample)
    data_plus_one_tissue = data.loc[idx_plus_one_tissue]
    
    idx_no_sample = [x for x in data.index if x != tissue_sample]
    
    data_minus_one_tissue = data.loc[idx_no_sample]
    
    return data_plus_one_tissue, data_minus_one_tissue

parser = argparse.ArgumentParser()
parser.add_argument('--gene', type=str, required=True)
parser.add_argument('--tissue', type=str, required=True)

if __name__ == '__main__':
    
    args = parser.parse_args()
    gene = args.gene
    tissue = args.tissue
    
    samples = pd.read_csv('../data/sample.tsv', sep='\t', index_col=0)
    participant_list = ['-'.join(x.split('-')[:2]) for x in samples.index]
    participants = pd.read_csv('../data/participant.tsv', sep='\t', index_col=0).loc[participant_list]
    annotation = pd.merge(samples, participants[['age', 'sex']], left_on='participant', right_index=True).drop_duplicates()

    data = pd.read_csv('Counts/' + gene + '.Counts.csv.gz', index_col=0)
    
    data_without_tissue = remove_tissue(data, annotation, tissue)
    data_without_tissue.to_csv('Counts/{gene}.Counts.no_{tissue}.csv.gz'.format(gene=gene, tissue=tissue))
    
    data_plus_one_tissue, data_minus_one_tissue = one_sample_for_tissue(data, annotation, tissue)
    data_plus_one_tissue.to_csv('Counts/{gene}.Counts.one_of_{tissue}.csv.gz'.format(gene=gene, tissue=tissue))
    
    data_minus_one_tissue.to_csv('Counts/{gene}.Counts.minus_one_of_{tissue}.csv.gz'.format(gene=gene, tissue=tissue))
    
    
    
    
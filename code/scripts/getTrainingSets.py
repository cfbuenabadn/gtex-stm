import numpy as np
import pandas as pd
import argparse
from sklearn.model_selection import train_test_split

def SplitTrainingSets(data, annotation):
    np.random.seed(0)
    annotation_data = annotation.loc[data.index]
    tissue_list = annotation_data.tissue_id.unique()
    
    idx_train = pd.Index([])
    idx_test = pd.Index([])
    
    for tissue in tissue_list:
        idx = annotation_data.loc[annotation_data.tissue_id == tissue].index
        train, test = train_test_split(idx, test_size=0.2)
        idx_train = idx_train.union(train)
        idx_test = idx_test.union(test)
        
    data_train = data.loc[idx_train]
    data_test = data.loc[idx_test]
    
    return data_train, data_test

parser = argparse.ArgumentParser()
parser.add_argument('--gene', type=str, required=True)

if __name__ == '__main__':
    
    args = parser.parse_args()
    gene = args.gene
    
    samples = pd.read_csv('../data/sample.tsv', sep='\t', index_col=0)
    participant_list = ['-'.join(x.split('-')[:2]) for x in samples.index]
    participants = pd.read_csv('../data/participant.tsv', sep='\t', index_col=0).loc[participant_list]
    annotation = pd.merge(samples, participants[['age', 'sex']], left_on='participant', right_index=True).drop_duplicates()

    data = pd.read_csv('Counts/' + gene + '.Counts.csv.gz', index_col=0)
    
    data_train, data_test = SplitTrainingSets(data, annotation)
    data_train.to_csv('Tests/{gene}/Counts/Train.csv.gz'.format(gene=gene))
    data_test.to_csv('Tests/{gene}/Counts/Test.csv.gz'.format(gene=gene))
    
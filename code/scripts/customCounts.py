import numpy as np
import pandas as pd
import argparse
from sklearn.model_selection import train_test_split


def subset_data(data, annotation, tissue_list, tissue_list_one):
    
    data_annotation = annotation.loc[data.index]
    
    subset_idx = data_annotation.loc[data_annotation.tissue_id.isin(tissue_list)].index
    
    np.random.seed(0)
    
    idx_one = []
    for tissue in tissue_list_one:
        idx_tissue = data_annotation.loc[data_annotation.tissue_id == tissue].index
        tissue_sample = np.random.choice(idx_tissue)
        idx_one.append(tissue_sample)
        
    idx_one = pd.Index(idx_one)
        
    subset_idx = subset_idx.union(idx_one)
    
    data_subset = data.loc[subset_idx]
    
    return data_subset, idx_one


if __name__ == '__main__':
    
    samples = pd.read_csv('../data/sample.tsv', sep='\t', index_col=0)
    participant_list = ['-'.join(x.split('-')[:2]) for x in samples.index]
    participants = pd.read_csv('../data/participant.tsv', sep='\t', index_col=0).loc[participant_list]
    annotation = pd.merge(samples, participants[['age', 'sex']], left_on='participant', right_index=True).drop_duplicates()

    data = pd.read_csv('Counts/PKM.Counts.csv.gz', index_col=0).iloc[:,:10100]

    data.to_csv('Tests/PKM/Counts/SubsetRegion.csv.gz')
    
    data_train, data_test = train_test_split(data, test_size=0.2)
    data_train.to_csv('Tests/PKM/Counts/SubsetRegion_train.csv.gz')
    data_test.to_csv('Tests/PKM/Counts/SubsetRegion_test.csv.gz')

    data_subset_tissue, idx_one = subset_data(data, annotation, ['Whole_Blood', 'Pancreas', 'Spleen', 'Kidney_Cortex', 'Liver'], [])
    data_subset_tissue.to_csv('Tests/PKM/Counts/SubsetTissues.csv.gz')

    data_train, data_test = train_test_split(data_subset_tissue, test_size=0.2)
    data_train.to_csv('Tests/PKM/Counts/SubsetTissues_train.csv.gz')
    data_test.to_csv('Tests/PKM/Counts/SubsetTissues_test.csv.gz')

    data_one_of_tissue, idx_one = subset_data(data, annotation, ['Whole_Blood', 'Pancreas', 'Spleen', 'Kidney_Cortex', 'Liver'], 
                                 ['Muscle_Skeletal', 'Brain_Cortex', 'Heart_Atrial_Appendage'])

    data_one_of_tissue.to_csv('Tests/PKM/Counts/SubsetTissues_plus_one.csv.gz')

    data_train_plus_one_idx = data_train.index.union(idx_one)

    data_train_plus_one = data.loc[data_train_plus_one_idx]

    data_train_plus_one.to_csv('Tests/PKM/Counts/SubsetTissues_train_plus_one.csv.gz')
    
    
    
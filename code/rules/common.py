import pandas as pd
import numpy as np
import os

genes = []
with open('config/genes.txt', 'r') as fh:
    for line in fh:
        genes.append(line.rstrip())
    
import json

def get_tissue_samples(annotation, tissue_id):
        
    female_samples = annotation.loc[(annotation.tissue_id == tissue_id) & (annotation.sex == 'Female')].index.sort_values()
    male_samples = annotation.loc[(annotation.tissue_id == tissue_id) & (annotation.sex == 'Male')].index.sort_values()
    
    return female_samples, male_samples


def select_samples(file_manifest, female_samples, male_samples, n=50):
    f = open(file_manifest)
    
    
    data = json.load(f)

    train = []
    train_count = 0
    for s in female_samples:
        for i in range(len(data)):
            if s == data[i]['file_name'].split('.')[0]:
                train.append(data[i]['file_name'].split('.')[0])
                train_count += 1
        if train_count == n:
            break


#     train_count = 0
    for s in male_samples:
        for i in range(len(data)):
            if s == data[i]['file_name'].split('.')[0]:
                train.append(data[i]['file_name'].split('.')[0])
                train_count += 1
        if train_count == (2*n):
            break
            
    return train

        
samples = pd.read_csv('../data/sample.tsv', sep='\t', index_col=0)
participant_list = ['-'.join(x.split('-')[:2]) for x in samples.index]
participants = pd.read_csv('../data/participant.tsv', sep='\t', index_col=0).loc[participant_list]

annotation = pd.merge(samples, participants[['age', 'sex']], left_on='participant', right_index=True).drop_duplicates()

brain_female, brain_male = get_tissue_samples(annotation, 'Brain_Cortex')
muscle_female, muscle_male = get_tissue_samples(annotation, 'Muscle_Skeletal')

brain_female_three_samples = select_samples('../data/file-manifest.json', brain_female, brain_male, n=3)
brain_male_three_samples = select_samples('../data/file-manifest.json', brain_female, brain_male, n=3)


# Note: there was a mistake in making female and male samples. The bug makes both subsets to have both male and female... 
# will correct later. This will have to be generalized when I add more tissues.

brain_female_three_samples_only = pd.Index(brain_female_three_samples).intersection(brain_female)
brain_male_three_samples_only = pd.Index(brain_male_three_samples).intersection(brain_male)


brain_female = brain_female.difference(pd.Index(brain_female_three_samples))
brain_male = brain_male.difference(pd.Index(brain_male_three_samples))


muscle_female_three_samples = select_samples('../data/file-manifest.json', muscle_female, muscle_male, n=3)
muscle_male_three_samples = select_samples('../data/file-manifest.json', muscle_female, muscle_male, n=3)

muscle_female_three_samples_only = pd.Index(muscle_female_three_samples).intersection(muscle_female)
muscle_male_three_samples_only = pd.Index(muscle_male_three_samples).intersection(muscle_male)

muscle_female = muscle_female.difference(pd.Index(muscle_female_three_samples))
muscle_male = muscle_male.difference(pd.Index(muscle_male_three_samples))

female_samples = brain_female_three_samples + muscle_female_three_samples
male_samples = brain_male_three_samples + muscle_male_three_samples
    
female_samples_only = list(brain_female_three_samples_only) + list(muscle_female_three_samples_only)
male_samples_only = list(brain_male_three_samples_only) + list(muscle_male_three_samples_only)
    
brain_train = select_samples('../data/file-manifest.json', brain_female, brain_male, n=50)
muscle_train = select_samples('../data/file-manifest.json', muscle_female, muscle_male, n=50)


train_samples = brain_train + muscle_train


brain_test = select_samples('../data/file-manifest.json', brain_female.difference(pd.Index(brain_train)), 
                            brain_male.difference(pd.Index(brain_train)), n=10)
muscle_test = select_samples('../data/file-manifest.json', muscle_female.difference(pd.Index(muscle_train)), 
                            muscle_male.difference(pd.Index(muscle_train)), n=10)

test_samples = brain_test + muscle_test


gtex_samples = train_samples + test_samples

with open('../data/test1.txt', 'w') as fh:
    for sample in gtex_samples:
        fh.write(sample + '\n')


import numpy as np
import pandas as pd
import json


def get_tissue_samples(annotation, tissue_id):
        
    female_samples = annotation.loc[(annotation.tissue_id == tissue_id) & (annotation.sex == 'Female')].index.sort_values()
    male_samples = annotation.loc[(annotation.tissue_id == tissue_id) & (annotation.sex == 'Male')].index.sort_values()
    
    return female_samples, male_samples


def select_samples(file_manifest, female_samples, male_samples, n=50):
    f = open(file_manifest)
    
    
    data = json.load(f)

    selected_json = []
    train = []
    train_count = 0
    for s in female_samples:
        for i in range(len(data)):
            if s == data[i]['file_name'].split('.')[0]:
                train.append(data[i]['file_name'].split('.')[0])
                selected_json.append(data[i])
                train_count += 1
        if train_count == n:
            break


#     train_count = 0
    for s in male_samples:
        for i in range(len(data)):
            if s == data[i]['file_name'].split('.')[0]:
                train.append(data[i]['file_name'].split('.')[0])
                selected_json.append(data[i])
                train_count += 1
        if train_count == (2*n):
            break
            
    return train, selected_json


def write_manifest(selected_json, file_manifest_out):

    with open(file_manifest_out, 'w') as outfile:
        json.dump(selected_json, outfile)
        
    

        
if __name__ == '__main__':

    samples = pd.read_csv('../data/sample.tsv', sep='\t', index_col=0)
    participant_list = ['-'.join(x.split('-')[:2]) for x in samples.index]
    participants = pd.read_csv('../data/participant.tsv', sep='\t', index_col=0).loc[participant_list]

    annotation = pd.merge(samples, participants[['age', 'sex']], left_on='participant', right_index=True).drop_duplicates()
    
    brain_female, brain_male = get_tissue_samples(annotation, 'Brain_Cortex')
    muscle_female, muscle_male = get_tissue_samples(annotation, 'Muscle_Skeletal')
    
    brain_female_three_samples, brain_female_three_samples_manifest = select_samples('../data/file-manifest.json', brain_female, brain_male, n=3)
    brain_male_three_samples, brain_male_three_samples_manifest = select_samples('../data/file-manifest.json', brain_female, brain_male, n=3)
    
    brain_female = brain_female.difference(pd.Index(brain_female_three_samples))
    brain_male = brain_male.difference(pd.Index(brain_male_three_samples))
    
    
    
    
    muscle_female_three_samples, muscle_female_three_samples_manifest = select_samples('../data/file-manifest.json', muscle_female, muscle_male, n=3)
    muscle_male_three_samples, muscle_male_three_samples_manifest = select_samples('../data/file-manifest.json', muscle_female, muscle_male, n=3)
    
    muscle_female = muscle_female.difference(pd.Index(muscle_female_three_samples))
    muscle_male = muscle_male.difference(pd.Index(muscle_male_three_samples))
    
    
    female_manifest = brain_female_three_samples_manifest + muscle_female_three_samples_manifest
    male_manifest = brain_male_three_samples_manifest + muscle_male_three_samples_manifest

    brain_train, brain_train_manifest = select_samples('../data/file-manifest.json', brain_female, brain_male, n=50)
    muscle_train, muscle_train_manifest = select_samples('../data/file-manifest.json', muscle_female, muscle_male, n=50)


    train_samples = brain_train + muscle_train


    brain_test, brain_test_manifest = select_samples('../data/file-manifest.json', brain_female.difference(pd.Index(brain_train)), 
                                brain_male.difference(pd.Index(brain_train)), n=10)
    muscle_test, muscle_test_manifest = select_samples('../data/file-manifest.json', muscle_female.difference(pd.Index(muscle_train)), 
                                muscle_male.difference(pd.Index(muscle_train)), n=10)

    test_samples = brain_test + muscle_test


    gtex_samples = train_samples + test_samples
    
    gtex_samples_manifest = brain_train_manifest + muscle_train_manifest + brain_test_manifest + muscle_test_manifest

    
    write_manifest(gtex_samples_manifest, 'gtex-download/file_manifest/TrainTest_2tissues.json')
    write_manifest(female_manifest, 'gtex-download/file_manifest/ThreeSamplesFemale.json')
    write_manifest(male_manifest, 'gtex-download/file_manifest/ThreeSamplesMale.json')
    
    
    with open('../data/test2.txt', 'w') as fh:
        for sample in gtex_samples:
            fh.write(sample + '\n')

      
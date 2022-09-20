import pandas as pd
import numpy as np
import os

genes = []
with open('config/genes.txt', 'r') as fh:
    for line in fh:
        genes.append(line.rstrip())
        
        
gtex_samples = pd.read_csv('config/samples.tsv', sep='\t', index_col=0)
tissue_list = ['Adipose_Subcutaneous', 'Muscle_Skeletal', 'Arterial_Tibial', 'Breast_Mammary_Tissue',
               'Skin_Not_Sun_Exposed_Suprapubic', 'Brain_Cortex', 'Thyroid', 'Lung', 'Spleen', 'Pancreas',
               'Heart_Atrial_Appendage', 'Heart_Left_Ventricle', 'Brain_Cerebellum', 'Cells_Cultured_fibroblasts', 
               'Whole_Blood', 'Artery_Aorta', 'Brain_Hypothalamus', 'Brain_Hippocampus', 'Liver', 'Kidney_Cortex'
              ]


rMATS_files = ['{ASE}.MATS.JCEC.txt', '{ASE}.MATS.JC.txt', 'fromGTF.{ASE}.txt', 'fromGTF.novelJunction.{ASE}.txt',
                  'fromGTF.novelSpliceSite.{ASE}.txt', 'JCEC.raw.input.{ASE}.txt', 'JC.raw.input.{ASE}.txt']
splicing_types = ['A3SS', 'A5SS', 'MXE', 'RI', 'SE']

brain_cortex_samples = gtex_samples.loc[gtex_samples.tissue_id == 'Brain_Cortex'].index
muscle_skeletal_samples = gtex_samples.loc[gtex_samples.tissue_id == 'Muscle_Skeletal'].index
whole_blood_samples = gtex_samples.loc[gtex_samples.tissue_id == 'Whole_Blood'].index
liver_samples = gtex_samples.loc[gtex_samples.tissue_id == 'Liver'].index

    
rsem_output_list = ['summary.txt']
for rf in rMATS_files:
    for ASE in splicing_types:
        new_file = rf.format(ASE=ASE)
        rsem_output_list.append(new_file)

def GetTissueList(wildcards):
    IndID_list = gtex_samples.loc[gtex_samples.tissue_id == wildcards.Tissue].index
    return IndID_list

def GetTissueListStr(wildcards):
    IndID_list = '.'.join(gtex_samples.loc[gtex_samples.tissue_id == wildcards.Tissue].index)
    return IndID_list

def GetTissueListStr2(wildcards):
    IndID_list = '.'.join(gtex_samples.loc[gtex_samples.tissue_id == wildcards.Tissue_2].index)
    return IndID_list

def GetTissueGroupList(wildcards):
    IndID_list = gtex_samples.loc[(gtex_samples.tissue_id == wildcards.Tissue) & (gtex_samples.group == wildcards.Group)].index
    return IndID_list

def GetTissueGroupListStr(wildcards):
    IndID_list = '.'.join(gtex_samples.loc[(gtex_samples.tissue_id == wildcards.Tissue) & (gtex_samples.group == wildcards.Group)].index)
    return IndID_list

def GetTissueGroupListStr2(wildcards):
    IndID_list = '.'.join(gtex_samples.loc[(gtex_samples.tissue_id == wildcards.Tissue_2) & (gtex_samples.group == wildcards.Group)].index)
    return IndID_list

def GetTissueBAM(wildcards):
    bams = expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam", 
        IndID = gtex_samples.loc[gtex_samples.tissue_id == wildcards.Tissue].index
        )
    
    return bams


def GetTissueJunc(wildcards):
    junc = expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc", 
        IndID = gtex_samples.loc[gtex_samples.tissue_id == wildcards.Tissue].index
        )
    
    return junc

    
def GetTissueGroupBAM(wildcards):
    bams = expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam", 
        IndID = gtex_samples.loc[(gtex_samples.tissue_id == wildcards.Tissue) & (gtex_samples.group == wildcards.Group)].index
        )
    
    return bams

def GetTissueGroupJunc(wildcards):
    junc = expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc", 
        IndID = gtex_samples.loc[(gtex_samples.tissue_id == wildcards.Tissue) & (gtex_samples.group == wildcards.Group)].index
        )
    
    return junc
    
def GetTissueGroupJunc2(wildcards):
    junc = expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue_2}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc", 
        IndID = gtex_samples.loc[(gtex_samples.tissue_id == wildcards.Tissue_2) & (gtex_samples.group == wildcards.Group)].index
        )
    
    return junc
    
    
# import json

# def get_tissue_samples(annotation, tissue_id):
        
#     female_samples = annotation.loc[(annotation.tissue_id == tissue_id) & (annotation.sex == 'Female')].index.sort_values()
#     male_samples = annotation.loc[(annotation.tissue_id == tissue_id) & (annotation.sex == 'Male')].index.sort_values()
    
#     f = open('../data/file-manifest.json')       
#     data = json.load(f)

#     female_samples_out = filter_with_manifest(data, female_samples)
#     male_samples_out = filter_with_manifest(data, male_samples)
#     return female_samples_out, male_samples_out


# def filter_with_manifest(data, samples_list):
#     samples_out = []
#     for sample in samples_list:
#         for i in range(len(data)):
#             if sample == data[i]['file_name'].split('.')[0]:
#                 samples_out.append(data[i]['file_name'].split('.')[0])

#     samples_out = sorted(samples_out)
#     return samples_out


# def get_train_and_test_groups_tissue(female_samples, male_samples):    
    
#     female_test = female_samples[:3]
#     train = female_samples[3:53] + male_samples[:50]
#     test = female_samples[53:63] + male_samples[50:60]  

#     return female_test, train, test

# def get_all_samples(tissue_id, return_sets=False):
#     female_samples, male_samples = get_tissue_samples(annotation, tissue_id)
#     female_test, train, test = get_train_and_test_groups_tissue(female_samples, male_samples)
    
#     if return_sets:
#         return female_test, train, test
    
#     all_samples = female_test + train + test
#     return all_samples

# # def GetTissueSamples(wildcards):
# #     tissue_samples = get_all_samples(wildcards.Tissue, False)
# #     return tissue_samples

# def GetTissueSamplesJoin(wildcards):
#     tissue_samples = get_all_samples(wildcards.Tissue, False)
#     tissue_samples = '.'.join(tissue_samples)
#     return tissue_samples


# samples = pd.read_csv('../data/sample.tsv', sep='\t', index_col=0)
# participant_list = ['-'.join(x.split('-')[:2]) for x in samples.index]
# participants = pd.read_csv('../data/participant.tsv', sep='\t', index_col=0).loc[participant_list]

# annotation = pd.merge(samples, participants[['age', 'sex']], left_on='participant', right_index=True).drop_duplicates()

# gtex_samples = ''
# brain_cortex_samples = get_all_samples('Brain_Cortex', False)
# muscle_skeletal_samples = get_all_samples('Muscle_Skeletal', False)


    

# # brain_cortex_female, brain_cortex_male = get_tissue_samples(annotation, 'Brain_Cortex')
# # muscle_skeletal_female, muscle_skeletal_male = get_tissue_samples(annotation, 'Muscle_Skeletal')

# # brain_cortex_female_test, brain_cortex_train, brain_cortex_test = get_train_and_test_groups_tissue(brain_cortex_female, brain_cortex_male)

# # brain_cortex = brain_cortex_female_test.union(brain_cortex_train).union(brain_cortex_test)

# # muscle_skeletal_female_test, muscle_skeletal_train, muscle_skeletal_test = get_train_and_test_groups_tissue(muscle_skeletal_female, muscle_skeletal_male)

# # muscle_skeletal = muscle_skeletal_female_test.union(muscle_skeletal_train).union(muscle_skeletal_test)





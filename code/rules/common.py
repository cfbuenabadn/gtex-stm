import pandas as pd
import numpy as np
import os

def more_mem_after_first_attempt(wildcards, attempt):
    if int(attempt) == 1:
        return 12000
    else:
        return 24000

test_list = ['Brain_Frontal_Cortex_BA9_v_Muscle_Skeletal',
             'Brain_Frontal_Cortex_BA9_v_Brain_Putamen_basal_ganglia',
             'Liver_v_Whole_Blood',
             'Lung_v_Skin_Not_Sun_Exposed_Suprapubic',
             'BA9_1_v_BA9_2'
             ]
    
prueba_genes = ['ENSG00000172009', 'ENSG00000112081', 'ENSG00000067225', 'ENSG00000132879']

top_genes = ['ENSG00000196531', 'ENSG00000196535', 'ENSG00000166925',
       'ENSG00000167642', 'ENSG00000125744', 'ENSG00000139116',
       'ENSG00000196586', 'ENSG00000170175', 'ENSG00000204681',
       'ENSG00000168710', 'ENSG00000139218', 'ENSG00000100030',
       'ENSG00000137094', 'ENSG00000110931', 'ENSG00000116754',
       'ENSG00000116857', 'ENSG00000174437', 'ENSG00000114796',
       'ENSG00000160469', 'ENSG00000120798', 'ENSG00000082258',
       'ENSG00000183020', 'ENSG00000135052', 'ENSG00000109184',
       'ENSG00000134250', 'ENSG00000151498', 'ENSG00000099814',
       'ENSG00000173575', 'ENSG00000072422', 'ENSG00000117625',
       'ENSG00000169020', 'ENSG00000157404', 'ENSG00000197321',
       'ENSG00000112531', 'ENSG00000165802', 'ENSG00000114098',
       'ENSG00000162735', 'ENSG00000114416', 'ENSG00000047188',
       'ENSG00000136295', 'ENSG00000124380', 'ENSG00000109917',
       'ENSG00000122359', 'ENSG00000125834', 'ENSG00000105048',
       'ENSG00000134313', 'ENSG00000077522', 'ENSG00000197386',
       'ENSG00000151729', 'ENSG00000168214']

genes_test = []
with open('config/genes.txt', 'r') as fh:
    for line in fh:
        genes_test.append(line.rstrip())
        
test_samples = list()
with open('../data/test_samples_selected.txt', 'r') as fh:
    for x in fh:
        test_samples.append(x.rstrip())


genes_negatives = []
with open('config/genes.Brain_Cortex_v_Muscle_Sleketal.negatives.txt', 'r') as fh:
    for line in fh:
        genes_negatives.append(line.rstrip())

genes_positives = []
with open('config/genes.Brain_Cortex_v_Muscle_Sleketal.positives.txt', 'r') as fh:
    for line in fh:
        genes_positives.append(line.rstrip())
        
genes_neg_pos = genes_negatives + genes_positives

quick_test_genes = ['NDUFA3', 'SRSF3', #'SRSF6', 
        'RTN2', 'SNRNP27', 'SPHK2', 'SUOX', 'TRAF3IP2', 'UCP3', #'UPF3A', 
        'RPL24', 'PRELID3B', 
        'PSAP', 'OAZ1', 'OAS1', 'RPS13', 'GPX3', 'HSP90AB1', "FBXO44"]

genes = sorted(set(genes_test + genes_neg_pos + quick_test_genes)) + ['ABI2']


# gtex_samples = pd.read_csv('config/samples.tsv', sep='\t', index_col=0)
gtex_samples = pd.read_csv('config/samples_extended.tsv', sep='\t', index_col=0)
tissue_list = ['Adipose_Subcutaneous', 'Muscle_Skeletal', 'Arterial_Tibial', 'Breast_Mammary_Tissue',
               'Skin_Not_Sun_Exposed_Suprapubic', 'Brain_Cortex', 'Thyroid', 'Lung', 'Spleen', 'Pancreas',
               'Heart_Atrial_Appendage', 'Heart_Left_Ventricle', 'Brain_Cerebellum', 'Cells_Cultured_fibroblasts', 
               'Whole_Blood', 'Artery_Aorta', 'Brain_Hypothalamus', 'Brain_Hippocampus', 'Liver', 'Kidney_Cortex',
               'Cells_EBV-transformed_lymphocytes', 'Brain_Frontal_Cortex_BA9', 'Brain_Anterior_cingulate_cortex_BA24',
               'Brain_Putamen_basal_ganglia', 'Brain_Caudate_basal_ganglia', 'Brain_Cerebellar_Hemisphere'
              ]


rMATS_files = ['{ASE}.MATS.JCEC.txt', '{ASE}.MATS.JC.txt', 'fromGTF.{ASE}.txt', 'fromGTF.novelJunction.{ASE}.txt',
                  'fromGTF.novelSpliceSite.{ASE}.txt', 'JCEC.raw.input.{ASE}.txt', 'JC.raw.input.{ASE}.txt']
splicing_types = ['A3SS', 'A5SS', 'MXE', 'RI', 'SE']

brain_cortex_samples = gtex_samples.loc[gtex_samples.tissue_id == 'Brain_Cortex'].index
muscle_skeletal_samples = gtex_samples.loc[gtex_samples.tissue_id == 'Muscle_Skeletal'].index
whole_blood_samples = gtex_samples.loc[gtex_samples.tissue_id == 'Whole_Blood'].index
liver_samples = gtex_samples.loc[gtex_samples.tissue_id == 'Liver'].index
brain_hippocampus_samples = gtex_samples.loc[gtex_samples.tissue_id == 'Brain_Hippocampus'].index
brain_hypothalamus_samples = gtex_samples.loc[gtex_samples.tissue_id == 'Brain_Hypothalamus'].index
brain_cerebellum_samples = gtex_samples.loc[gtex_samples.tissue_id == 'Brain_Cerebellum'].index

kidney_cortex_samples = gtex_samples.loc[gtex_samples.tissue_id == 'Kidney_Cortex'].index
lung_samples = gtex_samples.loc[gtex_samples.tissue_id == 'Lung'].index
heart_atrial_appendage_samples = gtex_samples.loc[gtex_samples.tissue_id == 'Heart_Atrial_Appendage'].index
spleen_samples = gtex_samples.loc[gtex_samples.tissue_id == 'Spleen'].index
skin_not_sun_exposed_suprapubic_samples = gtex_samples.loc[gtex_samples.tissue_id == 'Skin_Not_Sun_Exposed_Suprapubic'].index

LCL_samples = gtex_samples.loc[gtex_samples.tissue_id == 'Cells_EBV-transformed_lymphocytes'].index
fibroblast_samples = gtex_samples.loc[gtex_samples.tissue_id == 'Cells_Cultured_fibroblasts'].index

BA9_samples = gtex_samples.loc[gtex_samples.tissue_id == 'Brain_Frontal_Cortex_BA9'].index
BA24_samples =  gtex_samples.loc[gtex_samples.tissue_id == 'Brain_Anterior_cingulate_cortex_BA24'].index
putamen_samples = gtex_samples.loc[gtex_samples.tissue_id == 'Brain_Putamen_basal_ganglia'].index 
caudate_samples = gtex_samples.loc[gtex_samples.tissue_id == 'Brain_Caudate_basal_ganglia'].index 
cerebellarh_samples = gtex_samples.loc[gtex_samples.tissue_id == 'Brain_Cerebellar_Hemisphere'].index 


adipose_subcutaneous_samples = gtex_samples.loc[gtex_samples.tissue_id == 'Adipose_Subcutaneous'].index
Breast_Mammary_Tissue_samples = gtex_samples.loc[gtex_samples.tissue_id == 'Breast_Mammary_Tissue'].index 
Heart_Left_Ventricle_samples = gtex_samples.loc[gtex_samples.tissue_id == 'Heart_Left_Ventricle'].index  
Skin_Sun_Exposed_Lower_leg_samples = gtex_samples.loc[gtex_samples.tissue_id == 'Skin_Sun_Exposed_Lower_leg'].index  
Testis_samples = gtex_samples.loc[gtex_samples.tissue_id == 'Testis'].index  

    
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




# def GetTissueListStr_subsamples(wildcards):
#     IndID_list = '.'.join(gtex_samples.loc[gtex_samples.tissue_id == wildcards.Tissue].index)
#     return IndID_list

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

def GetTissueBai(wildcards):
    bams = expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.bai", 
        IndID = gtex_samples.loc[gtex_samples.tissue_id == wildcards.Tissue].index
        )
    
    return bams


def GetTissueJunc(wildcards):
    junc = expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc", 
        IndID = gtex_samples.loc[gtex_samples.tissue_id == wildcards.Tissue].index
        )
    
    return junc

def GetTissueJunc2(wildcards):
    junc = expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue_2}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc", 
        IndID = gtex_samples.loc[gtex_samples.tissue_id == wildcards.Tissue_2].index
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
    
tissue_list = ['Brain_Anterior_cingulate_cortex_BA24',
'Brain_Frontal_Cortex_BA9',  
'Heart_Atrial_Appendage',
'Lung',
'Skin_Not_Sun_Exposed_Suprapubic',
'Brain_Cortex',                          
'Brain_Putamen_basal_ganglia',  
'Liver',                   
'Muscle_Skeletal',  
'Whole_Blood']

new_tissues = ['Adipose_Subcutaneous', 'Breast_Mammary_Tissue', 'Heart_Left_Ventricle', 'Skin_Sun_Exposed_Lower_leg', 'Testis']

selected_genes = pd.read_csv('../data/protein_coding_genes.bed.gz', sep='\t', #'../data/selected_genes.bed', sep='\t', 
                             names = ['chrom', 'start', 'end', 'gene', 'gene_symbol', 'strand'])


genes_filtered = [x.rstrip() for x in open('../data/genes_filtered.txt', 'r').readlines()]

gregor_samples = pd.read_csv('config/gregor_samples_sberg.tsv', sep='\t', index_col=0)
gregor_samples = list(gregor_samples.index)
    
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

with open('../data/whole_genes.txt', 'r') as fh:
    whole_genes = [x.rstrip() for x in fh.readlines()]

tissue_sub_list = sorted(['Brain_Anterior_cingulate_cortex_BA24',
                      'Brain_Cortex',
                      'Brain_Frontal_Cortex_BA9',
                      'Brain_Putamen_basal_ganglia',
                      'Skin_Not_Sun_Exposed_Suprapubic',
                      'Liver',
                      'Lung', 
                      'Heart_Atrial_Appendage', 
                      'Muscle_Skeletal',
                      'Whole_Blood'])


subsamples = []
for tissue_ in tissue_sub_list:
    subsamples += list(gtex_samples.loc[(gtex_samples.tissue_id == tissue_) & (gtex_samples.group == 'train')].index[[0,25, 50, 75, 99]])
    
# for tissue_ in ['Brain_Frontal_Cortex_BA9', 'Muscle_Skeletal', 'Whole_Blood', 'Skin_Not_Sun_Exposed_Suprapubic']:
#     qqnorm_output.append(f'QTLs/{tissue_}/single_tissue.snmf_3.qqnorm.bed.gz')

def much_more_mem_after_first_attempt(wildcards, attempt):
    if int(attempt) == 1:
        return 4000
    else:
        return 52000

NumPvalsForPi1Chunks = 40



#### GWAS

colocs_df = pd.read_csv("config/ColocRunWildcards.tsv", index_col=0, comment='#', sep='\t', keep_default_na=False)
colocs_gwas = colocs_df.loc[colocs_df['FeatureCoordinatesRedefinedFor']=='ForGWASColoc']

def GetMolPhenotypesToColoc(wildcards):
    ProvidedMolPhenotypeList = colocs_df.loc[wildcards.ColocName]['MolPhenotypesToColoc']
    if ProvidedMolPhenotypeList == '':
        return ' '.join(PhenotypesToColoc)
    else:
        return ProvidedMolPhenotypeList


def GetColocTsvFormattedString(string):
    """
    return a snakemake input function that returns formatted string with wildcard values that match colocs_df based on ColocName index wildcard
    """
    def F(wildcards):
        return string.format(**colocs_df.loc[wildcards.ColocName].to_dict())
    return F



########### Gao data

gao_files = os.listdir('gao_data/bed_files/')

first_read_plus = pd.Index([x.split('.plus.')[0] for x in gao_files if x.endswith('.plus.bed.gz')])
first_read_minus = pd.Index([x.split('.minus.')[0] for x in gao_files if x.endswith('.minus.bed.gz')])

second_read_plus = pd.Index([x.split('.plus.second_read.')[0] for x in gao_files if x.endswith('.plus.second_read.bed.gz')])
second_read_minus = pd.Index([x.split('.minus.second_read.')[0] for x in gao_files if x.endswith('.minus.second_read.bed.gz')])

gao_complete_samples = sorted(first_read_plus.intersection(second_read_plus).intersection(first_read_minus).intersection(second_read_minus))

gao_genes = list(pd.read_csv('Annotations/gencode.v44.primary_assembly.protein_coding.genes.bed.gz', sep='\t').gene_id)

gwas_annot = pd.read_csv('config/gwas_annot_map.csv')

gwas_traits = [x.rstrip() for x in open('config/gwas_features.txt').readlines()]

gwas_traits_bjf79 = pd.read_csv("config/bjf79_gwas_table.tsv", index_col='gwas', sep='\t', 
                                usecols=['gwas', 'trait', 'FTPPath', "SummaryStatsLocalFilepath", 'ProcessingMethod', "Continuous"])



gwas_traits_for_coloc = [x for x in list(gwas_annot.trait) if x in gwas_traits]
gwas_traits_for_coloc += list(gwas_traits_bjf79.index)

# print("gwas_traits_for_coloc")
# print(len(gwas_traits_for_coloc))

#gwas_traits = [x for x in gwas_traits if x in list(gwas_annot.trait)]
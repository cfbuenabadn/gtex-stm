import numpy as np
import pandas as pd
from sklearn.preprocessing import scale
from scipy.stats import rankdata
from scipy.stats import norm
from sklearn.impute import SimpleImputer

def qqnorm(x):
    n=len(x)
    a=3.0/8.0 if n<=10 else 0.5
    return(norm.ppf( (rankdata(x)-a)/(n+1.0-2.0*a) ))

def get_PSI_table(bed, participant_list):
    bed_PSI = bed[bed.columns.intersection(participant_list)]
    bed_PSI.index = bed.pid
    bed_PSI = bed_PSI.replace('NA', np.nan).astype(float)
    bed_PSI = bed_PSI.loc[bed_PSI.isna().mean(axis=1) < 0.1]
    return bed_PSI

def impute_PSI(bed_PSI):
    imputer = SimpleImputer()
    imputed_df = pd.DataFrame(imputer.fit_transform(bed_PSI.T)).T
    imputed_df.columns = bed_PSI.columns 
    imputed_df.index = bed_PSI.index
    return imputed_df

def scale_PSI(imputed_df):
    df_scale = pd.DataFrame(scale(imputed_df, axis=1))
    df_scale.columns = imputed_df.columns
    df_scale.index = imputed_df.index
    return df_scale

def qqnorm_PSI(df_scale):

    df_qqnorm_list = []
    
    for sample in df_scale.columns:
        qqnorm_row = qqnorm(df_scale[sample])
        df_qqnorm_list.append(qqnorm_row)

    df_qqnorm_list = np.array(df_qqnorm_list).T
    df_qqnorm = pd.DataFrame(df_qqnorm_list)
    df_qqnorm.columns = df_scale.columns
    df_qqnorm.index = df_scale.index

    return df_qqnorm


def get_qqnorm_bed(bed, participant_list):
    
    bed6 = bed[bed.columns[:6]]
    bed_PSI = get_PSI_table(bed, participant_list)
    imputed_df = impute_PSI(bed_PSI)
    df_scale = scale_PSI(imputed_df)
    df_qqnorm = qqnorm_PSI(df_scale)
    qqnorm_bed = pd.merge(bed6, df_qqnorm, left_on = 'pid', right_index=True)
    
    return qqnorm_bed

def qqnorm_single_tissue(input_tissue, output_tissue, participant_list, K=3):
    ba9 = pd.read_csv(f'ebpmf_models/single_tissue/tables/{input_tissue}.EL.bed.gz', sep='\t')
    ba9_bed = get_qqnorm_bed(ba9, participant_list)
    return ba9_bed

def qqnorm_multi_tissue(snmf_bed, gtex_samples, tissue_id):
    tissue_samples = gtex_samples.loc[gtex_samples.tissue_id == tissue_id].index.intersection(snmf_bed.columns)
    all_samples = gtex_samples.index.intersection(snmf_bed.columns)
    multi_qqnorm = get_qqnorm_bed(snmf_bed, all_samples)
    bed_cols = list(multi_qqnorm.columns[:6])
    multi_qqnorm = multi_qqnorm[bed_cols + list(tissue_samples)]
    multi_qqnorm_cols = bed_cols + ['-'.join(x.split('-')[:2]) for x in multi_qqnorm.columns[6:]]
    multi_qqnorm.columns = multi_qqnorm_cols
    return multi_qqnorm


def wrapper_single_tissue(input_tissue_list, tissue_id_list, K = 3):
    participant_list = pd.read_csv('../data/participant.tsv', sep='\t', index_col=0).index
    for input_tissue, tissue_id in zip(input_tissue_list, tissue_id_list):
        print(tissue_id)
        tissue_qqnorm = qqnorm_single_tissue(input_tissue, tissue_id, participant_list, K=K)
        tissue_qqnorm.to_csv(f'QTLs/{tissue_id}/single_tissue.snmf_{str(K)}.qqnorm.bed.gz', sep='\t', index=False)

def wrapper_multi_tissue(tissue_list, K_list = [2, 3, 4, 5, 10]):
    gtex_samples = pd.read_csv('config/samples.tsv', sep='\t', index_col=0)
    for k in K_list:
        snmf_bed = pd.read_csv(f'ebpmf_models/filtered/snmf_{str(k)}/tables/EL.bed.gz', sep='\t')

        for tissue_id in tissue_list:
            print(tissue_id)
            tissue_qqnorm = qqnorm_multi_tissue(snmf_bed, gtex_samples, tissue_id)
            tissue_qqnorm.to_csv(f'QTLs/{tissue_id}/multi_tissue.snmf_{str(k)}.qqnorm.bed.gz', sep='\t', index=False)

if __name__ == '__main__':

    print('working on single tissue models')
    wrapper_single_tissue(['BA9', 'MS', 'WB', 'Skin'], 
                          ['Brain_Frontal_Cortex_BA9', 'Muscle_Skeletal', 'Whole_Blood', 'Skin_Not_Sun_Exposed_Suprapubic'])

    tissues = sorted(['Brain_Anterior_cingulate_cortex_BA24',
                      'Brain_Cortex',
                      'Brain_Frontal_Cortex_BA9',
                      'Brain_Putamen_basal_ganglia',
                      'Skin_Not_Sun_Exposed_Suprapubic',
                      'Liver',
                      'Lung', 
                      'Heart_Atrial_Appendage', 
                      'Muscle_Skeletal',
                      'Whole_Blood'])


    print('working on multi tissue models')
    wrapper_multi_tissue(tissues)



    
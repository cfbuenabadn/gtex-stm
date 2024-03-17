import pandas as pd
import numpy as np

if __name__ == '__main__':
    gtex_samples = pd.read_csv('config/samples.tsv', sep='\t', index_col=0)
    
    bam_list = list()
    with open('../data/test_bams.txt', 'r') as fh:
        for x in fh:
            bam_list.append(x.rstrip())

    tissue_samples = ['Brain_Putamen_basal_ganglia', 'Muscle_Skeletal', 'Liver', 'Whole_Blood',
                  'Skin_Not_Sun_Exposed_Suprapubic', 'Lung', 'Brain_Frontal_Cortex_BA9']
    gtex_samples = gtex_samples.loc[bam_list]


    rmats_samples = gtex_samples.loc[
    gtex_samples.tissue_id.isin(tissue_samples) & (
        (gtex_samples.sex=='female') | (gtex_samples.tissue_id=='Brain_Frontal_Cortex_BA9'))]

    junc_template = "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/TestSamples/juncs/"
    junc_template += "{IndID}.leafcutter.junc"

    def make_b(b_list):
        b = [junc_template.format(IndID = x) for x in b_list]
        return b

    paired_tests = [
        'Brain_Frontal_Cortex_BA9_v_Muscle_Skeletal',
        'Brain_Frontal_Cortex_BA9_v_Brain_Putamen_basal_ganglia',
        'Liver_v_Whole_Blood',
        'Lung_v_Skin_Not_Sun_Exposed_Suprapubic',
        'BA9_1_v_BA9_2'
    ]

    tissue_dir = {}

    for idx, df in rmats_samples.groupby(['tissue_id', 'sex']):
        tissue_name = idx[0]
        if idx[1] == 'female':
            all_samples = make_b([x for x in df.index if x != 'GTEX-11DXY-0011-R10b-SM-DO12C'])
            three_samples = make_b(df.iloc[[0, 2, 4]].index)
            tissue_dir.update({tissue_name:{'3Samples':three_samples, 
                                       '5Samples':all_samples}})
                    
    tissue_name = 'BA9_1'
    three_samples = make_b(rmats_samples.iloc[[0, 2, 5]].index)
    four_samples = make_b(rmats_samples.iloc[[0, 2, 4, 8]].index)

    tissue_dir.update({tissue_name:{'3Samples':three_samples, 
                                    '5Samples':four_samples}})

    tissue_name = 'BA9_2'
    three_samples = make_b(rmats_samples.iloc[[1, 3, 7]].index)
    four_samples = make_b(rmats_samples.iloc[[1, 3, 5, 7]].index)

    tissue_dir.update({tissue_name:{'3Samples':three_samples, 
                                    '5Samples':four_samples}})
    
    print('success with samples file\n')
    
    def write_line(x, tissue):
        
        sample = x.split('/')[-1].split('.')[0]
        
        path_to_kallisto = f'DifferentialSplicing/kallisto/quant/{sample}/'
        
        x_line = f'{sample}\t{tissue}\t{path_to_kallisto}\n'
        return x_line
    
    def write_tissue_lines(tissue1, tissue2, tissue_dir, n, fh):       
        fh.write('sample\tcondition\tpath\n')
        tissue1_samples = [write_line(x, tissue1) for x in tissue_dir.get(tissue1).get(n + 'Samples')]
        tissue2_samples = [write_line(x, tissue2) for x in tissue_dir.get(tissue2).get(n + 'Samples')]
        lines_samples = tissue1_samples + tissue2_samples
        fh.writelines(lines_samples)

    file_template = 'DifferentialSplicing/kallisto/sleuth/{tissue1}_v_{tissue2}_{n}Samples/input/groups_file.txt'

    for paired_test in paired_tests:
        tissue1, tissue2 = paired_test.split('_v_')

        with open(file_template.format(tissue1=tissue1, tissue2=tissue2, n='3'), 'w') as fh:
            write_tissue_lines(tissue1, tissue2, tissue_dir, '3', fh)
            

        with open(file_template.format(tissue1=tissue1, tissue2=tissue2, n='5'), 'w') as fh:
            write_tissue_lines(tissue1, tissue2, tissue_dir, '5', fh)

            
        print('success with ' + paired_test)











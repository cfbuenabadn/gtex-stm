import numpy as np
import pandas as pd

if __name__ == '__main__':
    gtex_samples = pd.read_csv('config/samples.tsv', sep='\t', index_col=0)
    
    # bam_list = list()
    # with open('../data/test_bams.txt', 'r') as fh:
    #     for x in fh:
    #         bam_list.append(x.rstrip())

    tissue_samples = ['Brain_Frontal_Cortex_BA9', 'Muscle_Skeletal', 'Liver', 'Whole_Blood']
    
    # tissue_samples = ['Brain_Putamen_basal_ganglia', 'Muscle_Skeletal', 'Liver', 'Whole_Blood',
    #               'Skin_Not_Sun_Exposed_Suprapubic', 'Lung', 'Brain_Frontal_Cortex_BA9']
    # gtex_samples = gtex_samples.loc[bam_list]
    
            
    rmats_samples = gtex_samples.loc[
    gtex_samples.tissue_id.isin(tissue_samples) & (
        (gtex_samples.sex=='female') | (gtex_samples.tissue_id=='Brain_Frontal_Cortex_BA9'))]

    bam_template = "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download-for-DS/{Tissue}/bams/"#"/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/TestSamples/bams/"
    bam_template += "{IndID}.Aligned.sortedByCoord.out.patched.md.bam"

    def make_b(b_list, tissue, bam_template):
        b = [bam_template.format(Tissue = tissue, IndID = x) for x in b_list]
        b = ','.join(b)
        return b

    for idx, df in rmats_samples.groupby(['tissue_id', 'sex']):
        if idx[1] == 'female':
            three_b = list(df.iloc[[0, 2, 4]].index)
            three_b = make_b(three_b)

            five_b = list(df.index)
            five_b = make_b(five_b)

            with open(f'DifferentialSplicing/rMATS/bfiles/{idx[0]}.3b.txt', 'w') as fh:
                fh.write(three_b)

            with open(f'DifferentialSplicing/rMATS/bfiles/{idx[0]}.5b.txt', 'w') as fh:
                fh.write(five_b)

    with open(f'DifferentialSplicing/rMATS/bfiles/BA9_1.3b.txt', 'w') as fh:
        b_list = make_b(rmats_samples.iloc[[0, 2, 5]].index)
        fh.write(b_list)

    with open(f'DifferentialSplicing/rMATS/bfiles/BA9_2.3b.txt', 'w') as fh:
        b_list = make_b(rmats_samples.iloc[[1, 3, 7]].index)
        fh.write(b_list)

    with open(f'DifferentialSplicing/rMATS/bfiles/BA9_1.5b.txt', 'w') as fh:
        b_list = make_b(rmats_samples.iloc[[0, 2, 4, 8]].index)
        fh.write(b_list)

    with open(f'DifferentialSplicing/rMATS/bfiles/BA9_2.5b.txt', 'w') as fh:
        b_list = make_b(rmats_samples.iloc[[1, 3, 5, 7]].index)
        fh.write(b_list)


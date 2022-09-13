import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--test', type=str, required=True)
parser.add_argument('--ds_tool', type=str, required=True)
parser.add_argument('--condition1', type=str, required=True)
parser.add_argument('--condition2', type=str, required=True)
parser.add_argument('--condition1_list', type=str, required=True)
parser.add_argument('--condition2_list', type=str, required=True)

if __name__ == '__main__':
    args = parser.parse_args()
    test = args.test
    ds_tool = args.ds_tool
    condition1 = args.condition1
    condition2 = args.condition2
    condition1_list = args.condition1_list.split('.')
    condition2_list = args.condition2_list.split('.')
    
    bam_template = "/project2/yangili1/cfbuenabadn/gtex-stm/code/gtex-download/bams/{Test}/{{IndID}}.Aligned.sortedByCoord.out.patched.md.bam"
    bam_template = bam_template.format(Test = test)
    
    if ds_tool == 'rMATS':
        
        b1 = ','.join([bam_template.format(IndID=x) for x in condition1_list])
        b2 = ','.join([bam_template.format(IndID=x) for x in condition2_list])
        
        with open("DifferentialSplicing/rMATS/{Test}/b1.txt".format(Test=test), 'w') as fh:
            fh.write(b1)
            
        with open("DifferentialSplicing/rMATS/{Test}/b2.txt".format(Test=test), 'w') as fh:
            fh.write(b2)
            
    elif ds_tool == 'leafcutter':
        
        group_df = pd.DataFrame()
        
        all_samples = condition1_list + condition2_list
        
        samples_col = ['{IndID}.Aligned.sortedByCoord.out.patched.md.bam'.format(IndID=IndID) for IndID in all_samples]
        condition_col = ([condition1]*len(condition1_list)) + ([condition2]*len(condition2_list))
        
        group_df['samples'] = samples_col
        group_df['condition'] = condition_col
        
        group_df.to_csv('DifferentialSplicing/leafcutter/{Test}/input/groups_file.txt'.format(Test=test), sep='\t',
                        index=False, header=False)
            
            
            
            
            
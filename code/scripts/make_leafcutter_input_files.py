import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--group', type=str, required=True)
parser.add_argument('--tissue1', type=str, required=True)
parser.add_argument('--tissue2', type=str, required=True)
parser.add_argument('--tissue1_list', type=str, required=True)
parser.add_argument('--tissue2_list', type=str, required=True)
parser.add_argument('--output', type=str, required=True)

if __name__ == '__main__':
    args = parser.parse_args()
    group = args.group
    tissue1 = args.tissue1
    tissue1_list = args.tissue1_list.split('.')
    tissue2 = args.tissue2
    tissue2_list = args.tissue2_list.split('.')
    output = args.output
        
    group_df = pd.DataFrame()

    all_samples = tissue1_list + tissue2_list

    samples_col = ['{IndID}.Aligned.sortedByCoord.out.patched.md.bam'.format(IndID=IndID) for IndID in all_samples]
    condition_col = ([tissue1]*len(tissue1_list)) + ([tissue2]*len(tissue2_list))

    group_df['samples'] = samples_col
    group_df['condition'] = condition_col

    group_df.to_csv(output, sep='\t',
                    index=False, header=False)

            
            
            
            

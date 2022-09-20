import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--group', type=str, required=True)
parser.add_argument('--tissue', type=str, required=True)
parser.add_argument('--tissue_list', type=str, required=True)
parser.add_argument('--output', type=str, required=True)

if __name__ == '__main__':
    args = parser.parse_args()
    group = args.group
    tissue = args.tissue
    tissue_list = args.tissue_list.split('.')
    output = args.output
    
    bam_template = "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{Tissue}/bams/{{IndID}}.Aligned.sortedByCoord.out.patched.md.bam"
    bam_template = bam_template.format(Tissue=tissue)
    
        
    b = ','.join([bam_template.format(IndID=x) for x in tissue_list])

    with open(output, 'w') as fh:
        fh.write(b)
            
        
    
            
            
            

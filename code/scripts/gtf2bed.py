import numpy as np
import pandas as pd
from gtfparse import read_gtf
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--gtf', type=str, required=True)

parser.add_argument('--output', type=str, required=True)

parser.add_argument('--extension', type=int, required=True)

if __name__ == '__main__':
    
    args = parser.parse_args()
    gtf_file = args.gtf
    extension = int(args.extension)
    output = args.output
        
    gtf = read_gtf(gtf_file)
    
    gtf_select = gtf.loc[(gtf.gene_type == 'protein_coding')&(gtf.feature=='gene')&np.array([x[:3]=='chr' for x in gtf.seqname])]
    gtf_bed = gtf_select[['seqname', 'start', 'end', 'gene_name']].copy()
    gtf_bed.start += (-extension)
    gtf_bed.end += extension
    
    gtf_bed.to_csv(output, sep='\t', index=False, header=False)
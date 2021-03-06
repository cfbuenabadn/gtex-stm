import argparse
import gzip
import os
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--output', type=str, required=True)
#parser.add_argument('--gene_dir', type=str, required=True)
parser.add_argument('FILES', type=str, nargs="+")

if __name__ == '__main__':
    
#     if not os.path.isdir('Counts'):
#         os.mkdir('Counts')

    args = parser.parse_args()
    
    output = args.output
    #gene_dir = args.gene_dir
    
    #files = [x for x in sorted(os.listdir(gene_dir)) if 'sorted' in x]
    files = args.FILES

    with gzip.open(output, 'wb') as fh_out:
        
        processed_files = 0
        
        for bed in files:
            IndID = bed.split('/')[-1].split('.')[0]
            
            position = []
            coverage = []
            with gzip.open(bed, 'rt') as fh:
                
                for row in fh:
                    row = row.rstrip().split('\t')
                    row_counts = [row[7]]*int(row[8])
                    coverage.extend(row_counts)
                    
                    if row[0][:3] == 'chr':
                        chrom = row[0]
                    else:
                        chrom = 'chr' + row[0]
                    
                    start = int(row[1]) + 1
                    
            if processed_files == 0:
                positions = np.arange(len(coverage)) + start
                
                fh_out.write(('IndID,'+ ','.join([chrom + ':' + str(x) for x in positions]) +'\n').encode())
                
            token = np.any([int(x) > 0 for x in coverage])
                
            if token:
                fh_out.write((IndID + ',' + ','.join(coverage) + '\n').encode())
            
            processed_files += 1
                
                    
            


            

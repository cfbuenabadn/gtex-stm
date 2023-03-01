import pandas as pd
import sys
import re

_, gene = sys.argv
df = pd.read_csv('Annotations/genes.tab.gz', sep='\t', names = ['ensembl_id', 'coord'], index_col=0)
gene_coord = re.split('[: -]', df.loc[gene].coord)

fh = open("coverage/tmp/{gene}.bed".format(gene=gene), 'w')
fh.write('\t'.join(gene_coord) + '\t' + gene)
fh.close()


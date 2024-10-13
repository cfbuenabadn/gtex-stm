import pandas as pd
import numpy as np
from sklearn.impute import SimpleImputer
import sys
import gzip


arguments = sys.argv
snmf_exons_file = arguments[1]
EL_file = arguments[2]
output_file = arguments[3]

#snmf_exons_bed = '/project2/mstephens/cfbuenabadn/gtex-stm/code/ebpmf_models/filtered/snmf_10/tables/snmf.3prime_corrected.exons.sorted.bed.gz'
snmf_exons = pd.read_csv(snmf_exons_file, sep='\t', names = ['chrom', 'start', 'end', 'gene_id', 
                                                 'transcript_id', 'strand', 'factors', 'exon_id'])

transcripts = snmf_exons.groupby(['gene_id', 'transcript_id']).factors.first().reset_index()


EL = pd.read_csv(EL_file, sep='\t') #'ebpmf_models/filtered/snmf_10/tables/EL.bed.gz', sep='\t')
# samples = pd.read_csv('../data/sample.tsv', sep='\t')
samples_list = EL.columns[6:]
coords_cols = list(EL.columns[:6])

print(EL.shape)

samples_list = list(EL[samples_list].columns[EL[samples_list].isna().mean(axis=0) < 0.9])
print(coords_cols + samples_list)
EL = EL[coords_cols + samples_list]

# samples = samples.set_index('entity:sample_id')
# EL = EL.loc[EL.isna().sum(axis=1) < 100]


X = EL[samples_list]
print(X.shape)
imp = SimpleImputer(missing_values=np.nan, strategy='mean')

X = imp.fit_transform(X)

X = pd.DataFrame(X)
X.columns = samples_list
X.index = EL.pid
X['gid'] = list(EL.gid)


fh = gzip.open(output_file, 'wb') #'ebpmf_models/filtered/snmf_10/tables/transcript.EL.bed.gz', 'wb')

first_line = ('\t'.join(list(EL.columns)) + '\n').encode()
fh.write(first_line)


for idx, df in EL.groupby('gid'):

    idx = idx.split('.')[0]

    first_col = [str(x) for x in list(df[df.columns[:6]].iloc[0])]
    
    df = df[samples_list].T
    gene_factor_list = sorted(df.columns)
    # df['tissue_id'] = list(samples.loc[samples_list].tissue_id)

    transcript_slice = transcripts.loc[transcripts.gene_id == idx]

    if transcript_slice.shape[0] == 1:
        continue

    transcript_EL = pd.DataFrame()

    for _, transcript_row in transcript_slice.iterrows():
        transcript_id = transcript_row.transcript_id
        factors = transcript_row.factors.split(':')
        factors_idx = [idx + '.' + x for x in factors]
        if len(factors_idx) == 1:
            transcript_series = X.loc[factors_idx[0], samples_list]
        else:
            # break
            transcript_series = X.loc[factors_idx, samples_list].sum(axis=0)

        transcripts_first_cols = first_col[:3] + [transcript_id, idx, first_col[-1]]
        transcript_series = transcripts_first_cols + [str(round(x, 4)) for x in transcript_series]

        transcript_line = ('\t'.join(transcript_series) + '\n').encode()
        fh.write(transcript_line)

fh.close()

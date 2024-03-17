import pandas as pd
gregor_juncs = pd.read_csv('gregor_data/sberger/processed_data/junc.bed.gz', sep='\t', 
                           names = ['chrom', 'start', 'end', 'sample', 'counts'])

chrom_list = ['chr' + str(x) for x in range(1, 23)] 
chrom_list += ['X', 'Y']
gregor_juncs = gregor_juncs[gregor_juncs.chrom.isin(chrom_list)]

gregor_juncs['junc_name'] = gregor_juncs.chrom + ':' + gregor_juncs.start.astype(str) + ':' + gregor_juncs.end.astype(str)

gregor_juncs = gregor_juncs.groupby(['junc_name', 'sample']).counts.sum().reset_index().pivot(index='junc_name', columns='sample', values='counts').fillna(0)

samples = gregor_juncs.columns
chrom = [x.split(':')[0] for x in gregor_juncs.index]
start = [x.split(':')[1] for x in gregor_juncs.index]
end = [x.split(':')[2] for x in gregor_juncs.index]

gregor_juncs['#chrom'] = chrom
gregor_juncs['start'] = start
gregor_juncs['end'] = end

gregor_juncs = gregor_juncs[['#chrom', 'start', 'end'] + samples]

gregor_juncs.to_csv('gregor_data/sberger/processed_data/junc.bed.gz', sep='\t', index=False, header=True)
#gregor_juncs.loc[(gregor_juncs > 0).astype(int).sum(axis=1) > ]
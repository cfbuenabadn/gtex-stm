import numpy as np
import pandas as pd
from pybedtools import BedTool
import gzip

snmf_exons_file = '/project2/mstephens/cfbuenabadn/gtex-stm/code/ebpmf_models/filtered/snmf_10/tables/snmf.merged_isoforms.exons.sorted.bed.gz'
snmf_exons_df = pd.read_csv(snmf_exons_file, sep='\t', names = ['chrom', 'start', 'end', 'gene_id', 
                                                 'transcript_id', 'strand', 'factors', 'exon_id'])

gencode_exons_file = '/project2/mstephens/cfbuenabadn/gtex-stm/code/Annotations/gencode.v44.primary_assembly.exons.sorted.bed.gz'
gencode_exons_df = pd.read_csv(gencode_exons_file, sep='\t', 
                            names = ['chrom', 'start', 'end', 'gene_id', 
                                                           'transcript_id', 'strand', 'exon_id', 'transcript_support_level',
                                                           'basic', 'Ensembl_canonical', 'MANE_Select', 'appris', 'transcript_type'])

snmf_introns_df = pd.read_csv('snmf_introns.bed.gz', sep='\t')
gencode_introns_df = pd.read_csv('gencode_introns.bed.gz', sep='\t')

snmf_introns_bed = BedTool.from_dataframe(snmf_introns_df.loc[snmf_introns_df.end >= snmf_introns_df.start]).sort()
gencode_introns_bed = BedTool.from_dataframe(gencode_introns_df.loc[gencode_introns_df.end >= gencode_introns_df.start]).sort()


annot_columns = ['chrom', 'start', 'end', 'gene_id', 'transcript_id', 'strand']

gencode_exons_df = gencode_exons_df.loc[(gencode_exons_df.end >= gencode_exons_df.start), annot_columns].drop_duplicates()
gencode_exons_bed = BedTool.from_dataframe(gencode_exons_df).sort()

snmf_exons_df = snmf_exons_df.loc[(snmf_exons_df.end >= snmf_exons_df.start), annot_columns].drop_duplicates()
snmf_exons_bed = BedTool.from_dataframe(snmf_exons_df).sort()

vastdb_exons_df = pd.read_csv('/project2/yangili1/cfbuenabadn/CSQ_additional_files/EVENT_INFO-hg38.tab.gz', sep='\t')
vastdb_exons_df = vastdb_exons_df.dropna()
vastdb_exons_cols = list(vastdb_exons_df.columns)
vastdb_exons_df['chrom'] = vastdb_exons_df.CO_A.apply(lambda x: x.split(':')[0])
vastdb_exons_df['start'] = vastdb_exons_df.CO_A.apply(lambda x: x.split(':')[1].split('-')[0])
vastdb_exons_df['end'] = vastdb_exons_df.CO_A.apply(lambda x: x.split('-')[1])
vastdb_exons_cols = ['chrom', 'start', 'end'] + vastdb_exons_cols
vastdb_exons_df = vastdb_exons_df[vastdb_exons_cols]

vastdb_exons_df = vastdb_exons_df.loc[(~vastdb_exons_df.COMPLEX.isin(['MIC-S', 'MIC_S', 'IR', 'Alt3', 'Alt5', 'MIC-M'])) & (vastdb_exons_df.EVENT.apply(lambda x: x.startswith('HsaEX')))]
vastdb_exons_df = vastdb_exons_df.loc[(~ vastdb_exons_df.COMPLEX.apply(lambda x: (x.startswith('ME') or (x == 'A_NA') or (x == 'ANN') or x == 'NA*')))]
vastdb_exons_df = vastdb_exons_df[[x for x in vastdb_exons_df.columns if x not in ['Seq_C1', 'Seq_A', 'Seq_C2']]]
vastdb_exons_df['strand'] = vastdb_exons_df.REF_CO.apply(lambda x: x.split(':')[-1])

intron_start_list = []
intron_end_list = []
for idx, row in vastdb_exons_df.iterrows():
    strand = row.strand
    if strand == '+':
        intron_start = int(row.CO_C1.split('-')[-1])
        intron_end = int(row.CO_C2.split('-')[0].split(':')[1])
    else:
        intron_start = int(row.CO_C2.split('-')[-1])
        intron_end = int(row.CO_C1.split('-')[0].split(':')[1])
        
    intron_start_list.append(intron_start)
    intron_end_list.append(intron_end)

vastdb_introns_df = pd.DataFrame()
vastdb_introns_df['chrom'] = list(vastdb_exons_df.chrom)
vastdb_introns_df['start'] = intron_start_list
vastdb_introns_df['end'] = intron_end_list
vastdb_introns_df['GENE'] = list(vastdb_exons_df.GENE)
vastdb_introns_df['EVENT'] = list(vastdb_exons_df.EVENT)
vastdb_introns_df['strand'] = list(vastdb_exons_df.strand)

vastdb_introns_bed = BedTool.from_dataframe(vastdb_introns_df).sort()
vastdb_exons_bed = BedTool.from_dataframe(vastdb_exons_df).sort()

introns_columns = ['chrom', 'start', 'end', 'gene_id', 'transcript_id', 'strand']


gencode_vastdb_exon_overlaps = gencode_exons_bed.intersect(vastdb_exons_bed, f=1, F=1, wo=True).to_dataframe(
    names = annot_columns + ['vastdb_' + x for x in list(vastdb_exons_df.columns)] + ['overlap']
)

gencode_vastdb_intron_overlaps = gencode_introns_bed.intersect(vastdb_introns_bed, f=1, F=1, wo=True).to_dataframe(
    names = introns_columns + ['vastdb_' + x for x in list(vastdb_introns_df.columns)] + ['overlap']
)

gencode_vastdb_exons_captured = pd.Index(gencode_vastdb_exon_overlaps.vastdb_EVENT.unique())
gencode_vastdb_introns_captured = pd.Index(gencode_vastdb_intron_overlaps.vastdb_EVENT.unique())

snmf_vastdb_exon_overlaps = snmf_exons_bed.intersect(vastdb_exons_bed, f=1, F=1, wo=True).to_dataframe(
    names = annot_columns + ['vastdb_' + x for x in list(vastdb_exons_df.columns)] + ['overlap']
)

snmf_vastdb_intron_overlaps = snmf_introns_bed.intersect(vastdb_introns_bed, f=1, F=1, wo=True).to_dataframe(
    names = introns_columns + ['vastdb_' + x for x in list(vastdb_introns_df.columns)] + ['overlap']
)

snmf_vastdb_exons_captured = pd.Index(snmf_vastdb_exon_overlaps.vastdb_EVENT.unique())
snmf_vastdb_introns_captured = pd.Index(snmf_vastdb_intron_overlaps.vastdb_EVENT.unique())


# def get_internal_exons(exons_df):
#     internal_exon_chrom_list = []
#     internal_exon_start_list = []
#     internal_exon_end_list = []
#     internal_exon_gene_list = []
#     internal_exon_transcript_list = []
#     internal_exon_strand_list = []
#     internal_exon_intron_start_list = []
#     internal_exon_intron_end_list = []
#     internal_exon_flanking1_start_list = []
#     internal_exon_flanking2_end_list = []
    
#     for transcript, df in tqdm(exons_df.groupby('transcript_id'), leave=True, position=0):
#         iso_len = df.shape[0]
#         if iso_len >= 3:
#             internal_exon_chrom = df.chrom.iloc[0]
#             internal_exon_gene = df.gene_id.iloc[0]
#             internal_exon_transcript = transcript
#             internal_exon_strand = df.strand.iloc[0]
#             for i in range(1, iso_len-1):
                
#                 internal_exon_start = df.iloc[i].start
#                 internal_exon_end = df.iloc[i].end
        
        
#                 internal_exon_intron_start = df.iloc[i-1].end
#                 internal_exon_intron_end = df.iloc[i+1].start

#                 internal_exon_flanking1_start = df.iloc[i-1].start
#                 internal_exon_flanking2_end = df.iloc[i+1].end
            
    
#                 internal_exon_chrom_list.append(internal_exon_chrom)
#                 internal_exon_start_list.append(internal_exon_start)
#                 internal_exon_end_list.append(internal_exon_end)
#                 internal_exon_gene_list.append(internal_exon_gene)
#                 internal_exon_transcript_list.append(internal_exon_transcript)
#                 internal_exon_strand_list.append(internal_exon_strand)
#                 internal_exon_intron_start_list.append(internal_exon_intron_start)
#                 internal_exon_intron_end_list.append(internal_exon_intron_end)

#                 internal_exon_flanking1_start_list.append(internal_exon_flanking1_start)
#                 internal_exon_flanking2_end_list.append(internal_exon_flanking2_end)
    
#     internal_exon_df = pd.DataFrame()
#     internal_exon_df['chrom'] = internal_exon_chrom_list
#     internal_exon_df['start'] = internal_exon_start_list
#     internal_exon_df['end'] = internal_exon_end_list
#     internal_exon_df['gene'] = internal_exon_gene_list
#     internal_exon_df['transcript'] = internal_exon_transcript_list
#     internal_exon_df['strand'] = internal_exon_strand_list
#     internal_exon_df['intron_start'] = internal_exon_intron_start_list
#     internal_exon_df['intron_end'] = internal_exon_intron_end_list

#     internal_exon_df['flanking1_start'] = internal_exon_flanking1_start_list
#     internal_exon_df['flanking2_end'] = internal_exon_flanking2_end_list


#     exon_coordinates = internal_exon_df.chrom + ':' + internal_exon_df.start.astype(str) + '-' + internal_exon_df.end.astype(str)
#     internal_exon_df['exon_coordinates'] = exon_coordinates
#     exon_name = internal_exon_df.exon_coordinates + ':' + internal_exon_df.gene
#     internal_exon_df['exon_id'] = exon_name

#     return internal_exon_df

# snmf_internal_exons_df = get_internal_exons(snmf_exons_df)
# gencode_internal_exons_df = get_internal_exons(gencode_exons_df)

# snmf_internal_exons_df = snmf_internal_exons_df[[x for x in snmf_internal_exons_df.columns if x != 'transcript']].drop_duplicates().reset_index(drop=True)
# gencode_internal_exons_df = gencode_internal_exons_df[[x for x in gencode_internal_exons_df.columns if x != 'transcript']].drop_duplicates().reset_index(drop=True)

# snmf_internal_exons_df.to_csv('snmf_internal_exons.bed.gz', sep='\t', header=True, index=False)
# gencode_internal_exons_df.to_csv('gencode_internal_exons.bed.gz', sep='\t', header=True, index=False)

snmf_internal_exons_df = pd.read_csv('snmf_internal_exons_w_transcript.bed.gz', sep='\t')
gencode_internal_exons_df = pd.read_csv('gencode_internal_exons_w_transcript.bed.gz', sep='\t')

snmf_internal_exons_df['transcripts'] = snmf_internal_exons_df.groupby(['chrom', 'start', 'end', 'gene', 'exon_coordinates', 'strand', 'intron_start', 'intron_end']).transcripts.transform(lambda x: '|'.join(sorted(x)))
gencode_internal_exons_df['transcripts'] = gencode_internal_exons_df.groupby(['chrom', 'start', 'end', 'gene', 'exon_coordinates', 'strand', 'intron_start', 'intron_end']).transcripts.transform(lambda x: '|'.join(sorted(x)))

snmf_internal_exons_df_flanking1 = pd.DataFrame(snmf_internal_exons_df.groupby(['chrom', 'start', 'end', 'gene', 'exon_coordinates', 'strand', 'intron_start', 'intron_end', 'transcripts']).flanking1_start.max().astype(int))
snmf_internal_exons_df_flanking2 = pd.DataFrame(snmf_internal_exons_df.groupby(['chrom', 'start', 'end', 'gene', 'exon_coordinates', 'strand', 'intron_start', 'intron_end', 'transcripts']).flanking2_end.min().astype(int))
snmf_internal_exons_df = snmf_internal_exons_df_flanking1.merge(snmf_internal_exons_df_flanking2, left_index=True, right_index=True).reset_index()

del snmf_internal_exons_df_flanking1
del snmf_internal_exons_df_flanking2

gencode_internal_exons_df_flanking1 = pd.DataFrame(gencode_internal_exons_df.groupby(['chrom', 'start', 'end', 'gene', 'exon_coordinates', 'strand', 'intron_start', 'intron_end', 'transcripts']).flanking1_start.max().astype(int))
gencode_internal_exons_df_flanking2 = pd.DataFrame(gencode_internal_exons_df.groupby(['chrom', 'start', 'end', 'gene', 'exon_coordinates', 'strand', 'intron_start', 'intron_end', 'transcripts']).flanking2_end.min().astype(int))
gencode_internal_exons_df = gencode_internal_exons_df_flanking1.merge(gencode_internal_exons_df_flanking2, left_index=True, right_index=True).reset_index()

del gencode_internal_exons_df_flanking1
del gencode_internal_exons_df_flanking2

snmf_internal_exons_bed = BedTool.from_dataframe(snmf_internal_exons_df).sort()
snmf_col_names = list(snmf_internal_exons_df.columns)
del snmf_internal_exons_df

gencode_internal_exons_bed = BedTool.from_dataframe(gencode_internal_exons_df).sort()
gencode_col_names = list(gencode_internal_exons_df.columns)
del gencode_internal_exons_df

del snmf_exons_df
del gencode_exons_df

del snmf_introns_df
del gencode_introns_df

intron_cols = [f'{x}_intron' for x in ['chrom', 'start', 'end', 'gene_id', 'transcript_id', 'strand']]
snmf_col_names += intron_cols #[f'{x}_intron' for x in snmf_introns_bed.to_dataframe().columns]
snmf_col_names += ['overlaps']

gencode_col_names += intron_cols
gencode_col_names += ['overlaps']

snmf_cassette_exons = snmf_internal_exons_bed.intersect(snmf_introns_bed, f=1, wo=True).to_dataframe(names = snmf_col_names)

# snmf_cassette_exons = snmf_internal_exons_intron_intersect.loc[snmf_internal_exons_intron_intersect.overlaps > 0]
snmf_cassette_exons = snmf_cassette_exons.loc[
(snmf_cassette_exons.intron_start == snmf_cassette_exons.start_intron) & (snmf_cassette_exons.intron_end == snmf_cassette_exons.end_intron)
][['chrom', 'start', 'end', 'gene', 'exon_coordinates', 'strand', 'intron_start', 'intron_end', 'flanking1_start', 'flanking2_end', 'transcripts']].drop_duplicates()

gencode_cassette_exons = gencode_internal_exons_bed.intersect(gencode_introns_bed, f=1, wo=True).to_dataframe(
    names = gencode_col_names)

# gencode_cassette_exons = gencode_internal_exons_intron_intersect.loc[gencode_internal_exons_intron_intersect.overlaps > 0]
gencode_cassette_exons = gencode_cassette_exons.loc[
(gencode_cassette_exons.intron_start == gencode_cassette_exons.start_intron) & (gencode_cassette_exons.intron_end == gencode_cassette_exons.end_intron)
][['chrom', 'start', 'end', 'gene', 'exon_coordinates', 'strand', 'intron_start', 'intron_end', 'flanking1_start', 'flanking2_end', 'transcripts']].drop_duplicates()

gencode_exon_name = []
current_gene = ''
for idx, row in gencode_cassette_exons.iterrows():
    gene_ = row.gene
    if gene_ != current_gene:
        gene_counter = 1
        current_gene = gene_
    exon_name = current_gene + ':' + str(gene_counter)
    gencode_exon_name.append(exon_name)
    gene_counter += 1

gencode_cassette_exons['exon_name'] = gencode_exon_name

gencode_cassette_exons_bed = BedTool.from_dataframe(gencode_cassette_exons)


gencode_cassette_exons_columns = gencode_cassette_exons.columns
snmf_cassette_exons_columns = snmf_cassette_exons.columns

snmf_gencode_exons_captured = snmf_internal_exons_bed.intersect(gencode_cassette_exons_bed, f=1, F=1, wo=True).to_dataframe(
    names = snmf_col_names[:11] + [f'{x}_gencode' for x in gencode_cassette_exons_columns] + ['overlaps']
)

snmf_gencode_exons_captured = snmf_gencode_exons_captured.loc[(snmf_gencode_exons_captured.intron_start.astype(int) == snmf_gencode_exons_captured.intron_start_gencode.astype(int)) & (snmf_gencode_exons_captured.intron_end.astype(int) == snmf_gencode_exons_captured.intron_end_gencode.astype(int))]

snmf_gencode_cassette = pd.Index(snmf_cassette_exons.exon_coordinates).intersection(pd.Index(snmf_gencode_exons_captured.exon_coordinates))

gencode_cassette_introns_bed = BedTool.from_dataframe(gencode_cassette_exons[['chrom', 'intron_start', 'intron_end', 'exon_name', 'start', 
                                                                              'end', 'flanking1_start', 'flanking2_end', 'transcripts']])

snmf_gencode_introns_captured = snmf_introns_bed.intersect(gencode_cassette_introns_bed, f=1, F=1, wo=True).to_dataframe(
    names = ['chrom', 'start', 'end', 'gene', 'transcript_id', 'strand', 'chrom_gencode', 'intron_start_gencode', 
             'intron_end_gencode', 'exon_name_gencode', 'start_gencode', 'end_gencode', 'flanking1_start_gencode', 'flanking2_end_gencode', 'transcripts', 'overlaps']
)

snmf_gencode_introns_captured = snmf_gencode_introns_captured[['chrom', 'start', 'end', 'gene', 'strand', 'chrom_gencode', 'intron_start_gencode', 
             'intron_end_gencode', 'exon_name_gencode', 'start_gencode', 'end_gencode', 'flanking1_start_gencode', 'flanking2_end_gencode', 'transcripts', 'overlaps']].drop_duplicates()

# Captured only Gencode exon
snmf_gencode_exons_captured.loc[~snmf_gencode_exons_captured.exon_coordinates.isin(snmf_gencode_cassette)]

# Captured cassette exon
snmf_gencode_exons_captured.loc[snmf_gencode_exons_captured.exon_coordinates.isin(snmf_gencode_cassette)]

# Captured intron only
snmf_gencode_introns_captured.loc[~snmf_gencode_introns_captured.exon_name_gencode.isin(snmf_gencode_exons_captured.exon_name_gencode)]

with open('CoveragePlots/bed_files/cassette_exons/gencode/exon_only.exon.bed', 'w') as fh:
    with open('CoveragePlots/bed_files/cassette_exons/gencode/exon_only.flanking1.bed', 'w') as fh1:
        with open('CoveragePlots/bed_files/cassette_exons/gencode/exon_only.flanking2.bed', 'w') as fh2:
            for idx, row in snmf_gencode_exons_captured.loc[~snmf_gencode_exons_captured.exon_coordinates.isin(snmf_gencode_cassette)].iterrows():
                exon_line = '\t'.join([row.chrom, str(row.start), str(row.end), row.exon_name_gencode, row.transcripts]) + '\n'
                flanking1_line = '\t'.join([row.chrom, str(row.flanking1_start_gencode), str(row.intron_start_gencode), row.exon_name_gencode, row.transcripts]) + '\n'
                flanking2_line = '\t'.join([row.chrom, str(row.intron_end_gencode), str(row.flanking2_end_gencode), row.exon_name_gencode, row.transcripts]) + '\n'
                fh.write(exon_line)
                fh1.write(flanking1_line)
                fh2.write(flanking2_line)
    
with open('CoveragePlots/bed_files/cassette_exons/gencode/cassette_exon.exon.bed', 'w') as fh:
    with open('CoveragePlots/bed_files/cassette_exons/gencode/cassette_exon.flanking1.bed', 'w') as fh1:
        with open('CoveragePlots/bed_files/cassette_exons/gencode/cassette_exon.flanking2.bed', 'w') as fh2:
            for idx, row in snmf_gencode_exons_captured.loc[snmf_gencode_exons_captured.exon_coordinates.isin(snmf_gencode_cassette)].iterrows():
                exon_line = '\t'.join([row.chrom, str(row.start), str(row.end), row.exon_name_gencode, row.transcripts]) + '\n'
                flanking1_line = '\t'.join([row.chrom, str(row.flanking1_start_gencode), str(row.intron_start_gencode), row.exon_name_gencode, row.transcripts]) + '\n'
                flanking2_line = '\t'.join([row.chrom, str(row.intron_end_gencode), str(row.flanking2_end_gencode), row.exon_name_gencode, row.transcripts]) + '\n'
                fh.write(exon_line)
                fh1.write(flanking1_line)
                fh2.write(flanking2_line)


gtex_gaps_df = pd.read_csv('coverage/GTEx_overlapping_regions.bed.gz', sep='\t')
gtex_gaps_bed = BedTool.from_dataframe(gtex_gaps_df)

introns_captured_cols = snmf_gencode_introns_captured.columns
introns_only_bed = BedTool.from_dataframe(snmf_gencode_introns_captured.loc[~snmf_gencode_introns_captured.exon_name_gencode.isin(snmf_gencode_exons_captured.exon_name_gencode)])
introns_only_filtered = introns_only_bed.intersect(gtex_gaps_bed, v=True).intersect(snmf_exons_bed, v=True).to_dataframe(names = introns_captured_cols)
# introns_only_filtered.iterrows()

with open('CoveragePlots/bed_files/cassette_exons/gencode/intron_only.exon.bed', 'w') as fh:
    with open('CoveragePlots/bed_files/cassette_exons/gencode/intron_only.flanking1.bed', 'w') as fh1:
        with open('CoveragePlots/bed_files/cassette_exons/gencode/intron_only.flanking2.bed', 'w') as fh2:
            for idx, row in introns_only_filtered.iterrows(): #snmf_gencode_introns_captured.loc[~snmf_gencode_introns_captured.exon_name_gencode.isin(snmf_gencode_exons_captured.exon_name_gencode)].iterrows():
                exon_line = '\t'.join([row.chrom, str(row.start_gencode), str(row.end_gencode), row.exon_name_gencode, row.transcripts]) + '\n'
                flanking1_line = '\t'.join([row.chrom, str(row.flanking1_start_gencode), str(row.intron_start_gencode), row.exon_name_gencode, row.transcripts]) + '\n'
                flanking2_line = '\t'.join([row.chrom, str(row.intron_end_gencode), str(row.flanking2_end_gencode), row.exon_name_gencode, row.transcripts]) + '\n'
                fh.write(exon_line)
                fh1.write(flanking1_line)
                fh2.write(flanking2_line)


with open('CoveragePlots/bed_files/cassette_exons/vastdb/exon_only.exon.bed', 'w') as fh:
    with open('CoveragePlots/bed_files/cassette_exons/vastdb/exon_only.flanking1.bed', 'w') as fh1:
        with open('CoveragePlots/bed_files/cassette_exons/vastdb/exon_only.flanking2.bed', 'w') as fh2:
            for idx, row in vastdb_exons_df.loc[vastdb_exons_df.EVENT.isin(snmf_vastdb_exons_captured.difference(snmf_vastdb_introns_captured))].iterrows():
                exon_line = '\t'.join([row.chrom, str(row.start), str(row.end), row.EVENT]) + '\n'
                if strand == '+':
                    flanking1_start = row.CO_C1.split(':')[1].split('-')[0]
                    flanking1_end = row.CO_C1.split(':')[1].split('-')[1]
                    flanking2_start = row.CO_C2.split(':')[1].split('-')[0]
                    flanking2_end = row.CO_C2.split(':')[1].split('-')[1]
                if strand == '-':
                    flanking1_start = row.CO_C2.split(':')[1].split('-')[0]
                    flanking1_end = row.CO_C2.split(':')[1].split('-')[1]
                    flanking2_start = row.CO_C1.split(':')[1].split('-')[0]
                    flanking2_end = row.CO_C1.split(':')[1].split('-')[1]
                    
                flanking1_line = '\t'.join([row.chrom, flanking1_start, flanking1_end, row.EVENT]) + '\n'
                flanking2_line = '\t'.join([row.chrom, flanking2_start, flanking2_end, row.EVENT]) + '\n'
                
                fh.write(exon_line)
                fh1.write(flanking1_line)
                fh2.write(flanking2_line)

with open('CoveragePlots/bed_files/cassette_exons/vastdb/cassette_exon.exon.bed', 'w') as fh:
    with open('CoveragePlots/bed_files/cassette_exons/vastdb/cassette_exon.flanking1.bed', 'w') as fh1:
        with open('CoveragePlots/bed_files/cassette_exons/vastdb/cassette_exon.flanking2.bed', 'w') as fh2:
            for idx, row in vastdb_exons_df.loc[vastdb_exons_df.EVENT.isin(snmf_vastdb_exons_captured.intersection(snmf_vastdb_introns_captured))].iterrows():
                exon_line = '\t'.join([row.chrom, str(row.start), str(row.end), row.EVENT]) + '\n'
                if strand == '+':
                    flanking1_start = row.CO_C1.split(':')[1].split('-')[0]
                    flanking1_end = row.CO_C1.split(':')[1].split('-')[1]
                    flanking2_start = row.CO_C2.split(':')[1].split('-')[0]
                    flanking2_end = row.CO_C2.split(':')[1].split('-')[1]
                if strand == '-':
                    flanking1_start = row.CO_C2.split(':')[1].split('-')[0]
                    flanking1_end = row.CO_C2.split(':')[1].split('-')[1]
                    flanking2_start = row.CO_C1.split(':')[1].split('-')[0]
                    flanking2_end = row.CO_C1.split(':')[1].split('-')[1]
                    
                flanking1_line = '\t'.join([row.chrom, flanking1_start, flanking1_end, row.EVENT]) + '\n'
                flanking2_line = '\t'.join([row.chrom, flanking2_start, flanking2_end, row.EVENT]) + '\n'

                fh.write(exon_line)
                fh1.write(flanking1_line)
                fh2.write(flanking2_line)


introns_captured_cols = vastdb_exons_df.columns
introns_only_bed = BedTool.from_dataframe(vastdb_exons_df.loc[vastdb_exons_df.EVENT.isin(snmf_vastdb_introns_captured.difference(snmf_vastdb_exons_captured))])
introns_only_filtered = introns_only_bed.intersect(gtex_gaps_bed, v=True).intersect(snmf_exons_bed, v=True).to_dataframe(names = introns_captured_cols)
# introns_only_filtered.iterrows()

with open('CoveragePlots/bed_files/cassette_exons/vastdb/intron_only.exon.bed', 'w') as fh:
    with open('CoveragePlots/bed_files/cassette_exons/vastdb/intron_only.flanking1.bed', 'w') as fh1:
        with open('CoveragePlots/bed_files/cassette_exons/vastdb/intron_only.flanking2.bed', 'w') as fh2:
            for idx, row in introns_only_filtered.iterrows(): #vastdb_exons_df.loc[vastdb_exons_df.EVENT.isin(snmf_vastdb_introns_captured.difference(snmf_vastdb_exons_captured))].iterrows():
                exon_line = '\t'.join([row.chrom, str(row.start), str(row.end), row.EVENT]) + '\n'
                if strand == '+':
                    flanking1_start = row.CO_C1.split(':')[1].split('-')[0]
                    flanking1_end = row.CO_C1.split(':')[1].split('-')[1]
                    flanking2_start = row.CO_C2.split(':')[1].split('-')[0]
                    flanking2_end = row.CO_C2.split(':')[1].split('-')[1]
                if strand == '-':
                    flanking1_start = row.CO_C2.split(':')[1].split('-')[0]
                    flanking1_end = row.CO_C2.split(':')[1].split('-')[1]
                    flanking2_start = row.CO_C1.split(':')[1].split('-')[0]
                    flanking2_end = row.CO_C1.split(':')[1].split('-')[1]
                    
                flanking1_line = '\t'.join([row.chrom, flanking1_start, flanking1_end, row.EVENT]) + '\n'
                flanking2_line = '\t'.join([row.chrom, flanking2_start, flanking2_end, row.EVENT]) + '\n'

                fh.write(exon_line)
                fh1.write(flanking1_line)
                fh2.write(flanking2_line)



exon_name_list = []
current_gene = ''
for idx, row in snmf_cassette_exons.iterrows():
    gene_ = row.gene
    if gene_ != current_gene:
        gene_counter = 1
        current_gene = gene_
    gene_name = current_gene + ':' + str(gene_counter)
    exon_name_list.append(gene_name)



snmf_cassette_exons['exon_name'] = exon_name_list

snmf_cassette_exons.to_csv('ebpmf_models/filtered/snmf_10/tables/cassette_exons.bed.gz', sep='\t', header=True, index=False)

snmf_gencode_cassette_exons = snmf_cassette_exons.loc[snmf_cassette_exons.exon_coordinates.isin(snmf_gencode_exons_captured.exon_coordinates)]


snmf_cassette_exons_bed = BedTool.from_dataframe(snmf_cassette_exons)

snmf_vastdb_cassette_exons = snmf_cassette_exons_bed.intersect(vastdb_exons_bed, f=1, F=1, wo=True).to_dataframe(
    names = list(snmf_cassette_exons.columns) + [f'{x}_vastdb' for x in vastdb_exons_df.columns] + ['overlaps']
)

# intron_start_vastdb = []
# intron_end_vastdb = []
# for idx, row in snmf_vastdb_cassette_exons.iterrows():
#     if row.strand_vastdb == '+':
#         intron_start = row.CO_C1_vastdb
#         intron_end = row.CO_C2_vastdb
#     else:
#         intron_start = row.CO_C2_vastdb
#         intron_end = row.CO_C1_vastdb
#     intron_start = int(intron_start.split('-')[-1])
#     intron_end = int(intron_end.split(':')[-1].split('-')[0])

#     intron_start_vastdb.append(intron_start)
#     intron_end_vastdb.append(intron_end)

# snmf_vastdb_cassette_exons['intron_start_vastdb'] = intron_start_vastdb
# snmf_vastdb_cassette_exons['intron_end_vastdb'] = intron_end_vastdb

# snmf_vastdb_cassette_exons = snmf_vastdb_cassette_exons.loc[
# (snmf_vastdb_cassette_exons.intron_start == snmf_vastdb_cassette_exons.intron_start_vastdb) & (snmf_vastdb_cassette_exons.intron_end == snmf_vastdb_cassette_exons.intron_end_vastdb)
# ]

def confirm_junctions(cassette_exons, junctions_bed):
    
    I1 = cassette_exons[['chrom', 'intron_start', 'start', 'gene', 'exon_coordinates']].copy()
    I1.columns = ['chrom', 'start', 'end', 'gene', 'exon_coordinates']
    I1.start += 1
    I1.end -= 1
    
    I2 = cassette_exons[['chrom', 'end', 'intron_end', 'gene', 'exon_coordinates']].copy()
    I2.columns = ['chrom', 'start', 'end', 'gene', 'exon_coordinates']
    I2.start += 1
    I2.end -= 1
    
    SE = cassette_exons[['chrom', 'intron_start', 'intron_end', 'gene', 'exon_coordinates']].copy()
    SE.columns = ['chrom', 'start', 'end', 'gene', 'exon_coordinates']
    SE.start += 1
    SE.end -= 1
    
    I1_ = BedTool.from_dataframe(I1).intersect(junctions_bed, f=1, F=1, wao=True).to_dataframe(names = list(I1.columns) + [f'{x}_junc' for x in junctions.columns] + ['overlaps'])
    I2_ = BedTool.from_dataframe(I2).intersect(junctions_bed, f=1, F=1, wao=True).to_dataframe(names = list(I2.columns) + [f'{x}_junc' for x in junctions.columns] + ['overlaps'])
    SE_ = BedTool.from_dataframe(SE).intersect(junctions_bed, f=1, F=1, wao=True).to_dataframe(names = list(SE.columns) + [f'{x}_junc' for x in junctions.columns] + ['overlaps'])
    
    I1_confirmed = pd.Index(I1_.loc[I1_.overlaps > 0].exon_coordinates.unique())
    I2_confirmed = pd.Index(I2_.loc[I2_.overlaps > 0].exon_coordinates.unique())
    SE_confirmed = pd.Index(SE_.loc[SE_.overlaps > 0].exon_coordinates.unique())

    return I1_confirmed, I2_confirmed, SE_confirmed


junctions = pd.read_csv('gtex_tables/junctions.bed.gz', sep='\t')
junctions['gene'] = junctions.gene.apply(lambda x: x.split('.')[0])
junctions_bed = BedTool.from_dataframe(junctions)
snmf_I1, snmf_I2, snmf_SE = confirm_junctions(snmf_cassette_exons, junctions_bed)

snmf_junctions_cassette_exons = snmf_cassette_exons.loc[snmf_cassette_exons.exon_coordinates.isin(snmf_I1.intersection(snmf_I2))]

snmf_unconfirmed_cassette_exons = snmf_cassette_exons.loc[~snmf_cassette_exons.exon_coordinates.isin(
    pd.Index(snmf_gencode_cassette_exons.exon_coordinates).union(pd.Index(snmf_vastdb_cassette_exons.exon_coordinates)).union(
        pd.Index(snmf_junctions_cassette_exons.exon_coordinates))
)]

with open('CoveragePlots/bed_files/cassette_exons/snmf/gencode.exon.bed', 'w') as fh:
    with open('CoveragePlots/bed_files/cassette_exons/snmf/gencode.flanking1.bed', 'w') as fh1:
        with open('CoveragePlots/bed_files/cassette_exons/snmf/gencode.flanking2.bed', 'w') as fh2:
            for idx, row in snmf_gencode_cassette_exons.iterrows():
                exon_line = '\t'.join([row.chrom, str(row.start), str(row.end), row.exon_name, row.transcripts]) + '\n'
                flanking1_line = '\t'.join([row.chrom, str(row.flanking1_start), str(row.intron_start), row.exon_name, row.transcripts]) + '\n'
                flanking2_line = '\t'.join([row.chrom, str(row.intron_end), str(row.flanking2_end), row.exon_name, row.transcripts]) + '\n'
                fh.write(exon_line)
                fh1.write(flanking1_line)
                fh2.write(flanking2_line)
    
with open('CoveragePlots/bed_files/cassette_exons/snmf/vastdb.exon.bed', 'w') as fh:
    with open('CoveragePlots/bed_files/cassette_exons/snmf/vastdb.flanking1.bed', 'w') as fh1:
        with open('CoveragePlots/bed_files/cassette_exons/snmf/vastdb.flanking2.bed', 'w') as fh2:
            for idx, row in snmf_vastdb_cassette_exons.iterrows():
                exon_line = '\t'.join([row.chrom, str(row.start), str(row.end), row.exon_name, row.transcripts]) + '\n'
                flanking1_line = '\t'.join([row.chrom, str(row.flanking1_start), str(row.intron_start), row.exon_name, row.transcripts]) + '\n'
                flanking2_line = '\t'.join([row.chrom, str(row.intron_end), str(row.flanking2_end), row.exon_name, row.transcripts]) + '\n'
                fh.write(exon_line)
                fh1.write(flanking1_line)
                fh2.write(flanking2_line)

with open('CoveragePlots/bed_files/cassette_exons/snmf/junctions.exon.bed', 'w') as fh:
    with open('CoveragePlots/bed_files/cassette_exons/snmf/junctions.flanking1.bed', 'w') as fh1:
        with open('CoveragePlots/bed_files/cassette_exons/snmf/junctions.flanking2.bed', 'w') as fh2:
            for idx, row in snmf_junctions_cassette_exons.iterrows():
                exon_line = '\t'.join([row.chrom, str(row.start), str(row.end), row.exon_name, row.transcripts]) + '\n'
                flanking1_line = '\t'.join([row.chrom, str(row.flanking1_start), str(row.intron_start), row.exon_name, row.transcripts]) + '\n'
                flanking2_line = '\t'.join([row.chrom, str(row.intron_end), str(row.flanking2_end), row.exon_name, row.transcripts]) + '\n'
                fh.write(exon_line)
                fh1.write(flanking1_line)
                fh2.write(flanking2_line)

with open('CoveragePlots/bed_files/cassette_exons/snmf/unmatched.exon.bed', 'w') as fh:
    with open('CoveragePlots/bed_files/cassette_exons/snmf/unmatched.flanking1.bed', 'w') as fh1:
        with open('CoveragePlots/bed_files/cassette_exons/snmf/unmatched.flanking2.bed', 'w') as fh2:
            for idx, row in snmf_unconfirmed_cassette_exons.iterrows():
                exon_line = '\t'.join([row.chrom, str(row.start), str(row.end), row.exon_name, row.transcripts]) + '\n'
                flanking1_line = '\t'.join([row.chrom, str(row.flanking1_start), str(row.intron_start), row.exon_name, row.transcripts]) + '\n'
                flanking2_line = '\t'.join([row.chrom, str(row.intron_end), str(row.flanking2_end), row.exon_name, row.transcripts]) + '\n'
                fh.write(exon_line)
                fh1.write(flanking1_line)
                fh2.write(flanking2_line)

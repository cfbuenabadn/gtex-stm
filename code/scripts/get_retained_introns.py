import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import pybedtools
import gzip


def extract_retained_introns(protein_coding_bed, retained_intron_bed):
    X = pybedtools.BedTool.from_dataframe(protein_coding_bed)
    Y = pybedtools.BedTool.from_dataframe(retained_intron_bed)

    intersection = Y.intersect(X, c=True).to_dataframe()
    intersection.columns = ['chrom', 'start', 'end', 'overlaps']
    Z = pybedtools.BedTool.from_dataframe(intersection.loc[intersection.overlaps >= 2])

    if len(Z) >= 1:

        Z = Z.subtract(X).to_dataframe()
        Z.start -= 1
        Z.end += 1
        Z = Z[['chrom', 'start', 'end']]
        Z = pybedtools.BedTool.from_dataframe(Z).intersect(X, c=True).to_dataframe()
        Z = Z.loc[Z.name >= 2]
    
        if len(Z) >= 1:
            Z.start += 1
            Z.end -= 1

    else:
        Z = pd.DataFrame(columns = ['chrom', 'start', 'end'])

    return Z

def get_coverage_percentage(bed, start, end, hard_start, hard_end):

    if len(bed) == 1:
        x = np.max([((merged_intron_retention.start - start)/(end - start)).iloc[0], 0])
        y = np.min([((merged_intron_retention.end - start)/(end - start)).iloc[0], 1])
        covered = np.abs(y-x)
    else:
    
        X = ((bed.end - start)/(end - start)).iloc[:-1]
        Y = ((bed.start - start)/(end - start)).iloc[1:]
    
        covered = 0
    
        if (X.iloc[0] > 0) and (X.iloc[0] < 1):
            x = X.iloc[0]
            covered += x
    
        for i in range(len(X)-1):
            x = np.min([np.max([0, X.iloc[i+1]]), 1])
            y = np.max([np.min([1, Y.iloc[i]]), 0])
            
            covered += np.abs(y-x)
    
        if (Y.iloc[-1] < 1) and (Y.iloc[-1] > 0):
            y = Y.iloc[-1]
            covered += np.abs(1-y)

    return covered

gencode_exons_bed = '/project2/mstephens/cfbuenabadn/gtex-stm/code/Annotations/gencode.v44.primary_assembly.exons.sorted.bed.gz'
gencode_exons = pd.read_csv(gencode_exons_bed, sep='\t', 
                            names = ['chrom', 'start', 'end', 'gene_id', 
                                                           'transcript_id', 'strand', 'exon_id', 'transcript_support_level',
                                                           'basic', 'Ensembl_canonical', 'MANE_Select', 'appris', 'transcript_type'])

snmf_exons_bed = '/project2/mstephens/cfbuenabadn/gtex-stm/code/ebpmf_models/filtered/snmf_10/tables/snmf.merged_isoforms.exons.sorted.bed.gz'
snmf_exons = pd.read_csv(snmf_exons_bed, sep='\t', names = ['chrom', 'start', 'end', 'gene_id', 
                                                 'transcript_id', 'strand', 'factors', 'exon_id'])


annotation = pd.read_csv('/project2/mstephens/cfbuenabadn/gtex-stm/code/ebpmf_models/filtered/snmf_10/tables/annotated.snmf.merged_isoforms.tab.gz', sep='\t')
annotation = annotation.loc[annotation.appris_transcript_length != 'appris_transcript_length']
annotation['appris_transcript_length'] = annotation.appris_transcript_length.astype(int)
annotation['gene_id'] = [x.split('.')[0] for x in annotation.transcript]


with gzip.open('CoveragePlots/bed_files/introns/snmf.retained_introns.bed.gz', 'wb') as fh:

    for gene, df in snmf_exons.groupby('gene_id'):

        
        
        annot_slice = annotation.loc[annotation.gene_id == gene]
        annotated_isoforms = list(annot_slice.loc[annot_slice.appris_ref.apply(lambda x: ('intron.chain' in x) or ('ref.transcript' in x))].transcript)
    
        retained_introns = list(annot_slice.loc[annot_slice.appris_ref.apply(lambda x: 'retained.intron' in x)].transcript)
    
        bedtool = pybedtools.BedTool.from_dataframe(df.loc[df.transcript_id.isin(annotated_isoforms)])
        merged_bedtool = bedtool.merge().to_dataframe()
        if (len(annotated_isoforms) > 0) and (len(merged_bedtool) >= 2):
            
            start = merged_bedtool.end.min()
            end = merged_bedtool.start.max()
            hard_start = df.start.min()
            hard_end = df.end.max()
            
            if (len(retained_introns) > 0):
                
                merged_intron_retention = pybedtools.BedTool.from_dataframe(df.loc[df.transcript_id.isin(retained_introns)])
                merged_intron_retention = merged_intron_retention.merge().to_dataframe()
    
                if len(merged_intron_retention) == 1:
                    percent = get_coverage_percentage(merged_intron_retention, start, end, hard_start, hard_end)
                else:
                    percent = 1
    
                if (len(merged_intron_retention) >= 1):# or ((len(merged_intron_retention) == 1) and (percent < 0.5)):
    
                    # try:

                    strand = df.iloc[0].strand
                    transcripts = '|'.join(retained_introns)
                    

                    retained_i = extract_retained_introns(merged_bedtool, merged_intron_retention)

                    for idx, row in retained_i.iterrows():
                        line = ('\t'.join([row.chrom, str(row.start), str(row.end), gene, transcripts, strand]) + '\n')
                        # print(line)
                        line = line.encode()
                        fh.write(line)
                        
                    # except:
                    #     continue
with gzip.open('CoveragePlots/bed_files/introns/gencode.retained_introns.bed.gz', 'wb') as fh:

    for gene, df in gencode_exons.loc[gencode_exons.gene_id.isin(snmf_exons.gene_id.unique())].groupby('gene_id'):
        # df = gencode_exons.loc[(gencode_exons.gene_id == gene)]
        
        annotated_isoforms = list(df.loc[(df.transcript_type == 'protein_coding') & (~df.appris.isna())].transcript_id)
    
        retained_introns = list(df.loc[df.transcript_type == 'retained_intron'].transcript_id)
    
        bedtool = pybedtools.BedTool.from_dataframe(df.loc[df.transcript_id.isin(annotated_isoforms)])
        merged_bedtool = bedtool.merge().to_dataframe()
        if (len(annotated_isoforms) > 0) and (len(merged_bedtool) >= 2):
        
            start = merged_bedtool.end.min()
            end = merged_bedtool.start.max()
            hard_start = df.start.min()
            hard_end = df.end.max()
            
            if (len(retained_introns) > 0):
    
                merged_intron_retention = pybedtools.BedTool.from_dataframe(df.loc[df.transcript_id.isin(retained_introns)])
                merged_intron_retention = merged_intron_retention.merge().to_dataframe()
    
                if len(merged_intron_retention) >= 1:
    
    
                    try:

                        strand = df.iloc[0].strand
                        transcripts = '|'.join(retained_introns)
    
                        retained_i = extract_retained_introns(merged_bedtool, merged_intron_retention)

                        for idx, row in retained_i.iterrows():
                            line = ('\t'.join([row.chrom, str(row.start), str(row.end), gene, transcripts, strand]) + '\n').encode()
                            fh.write(line)
        
                    except:
                        continue


with gzip.open('CoveragePlots/bed_files/introns/gencode.appris_introns.bed.gz', 'wb') as fh:

    appris_transcripts = gencode_exons.loc[
    (gencode_exons.transcript_type == 'protein_coding') & (gencode_exons.appris == 'appris_principal_1') & gencode_exons.gene_id.isin(snmf_exons.gene_id.unique())
    ]

    for transcript, df in appris_transcripts.groupby('transcript_id'):

        if len(df) >= 2:
            chrom = df.iloc[0].chrom
            strand = df.iloc[0].strand
            gene = df.iloc[0].gene_id
            
            intron_pos = zip(df.end.iloc[:-1], df.start.iloc[1:])

            for intron_start, intron_end in intron_pos:
                line = ('\t'.join([chrom, str(intron_start), str(intron_end), gene, transcript, strand]) + '\n').encode()
                fh.write(line)
    
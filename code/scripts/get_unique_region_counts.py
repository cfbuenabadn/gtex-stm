import numpy as np
import pandas as pd

import gzip
from pybedtools import BedTool

import sys



def get_exon_coords(row):
    chrom = row.chrom
    start_coords = int(row.start)
    end_coords = int(row.end)
    exon_coords = [chrom + ':' + str(x) for x in range(start_coords, end_coords+1)]
    return exon_coords

def get_transcript_coords(transcript_df):
    transcript_coords = []
    for idx, row in transcript_df.iterrows():
        exon_coords = get_exon_coords(row)
        transcript_coords.extend(exon_coords)

    transcript_coords = pd.Index(transcript_coords)
    return transcript_coords

def get_meta(gene_df):
    meta_bed = BedTool.from_dataframe(gene_df).merge()
    return meta_bed

def get_difference_coords(gene1_bed, gene2_bed):
    sub_bed = gene1_bed.subtract(gene2_bed).to_dataframe()
    sub_coords = get_transcript_coords(sub_bed)
    return sub_coords

def get_counts_per_bp(coords, counts):
    coords = counts.columns.intersection(coords)
    cpb = counts[coords].sum(axis=1)/len(coords)
    return cpb

def process_gene(snmf, gencode, gene_id, min_len = 100):
    bed_cols = ['chrom', 'start', 'end']
    snmf_ = snmf.loc[snmf.gene_id == gene_id, bed_cols]
    gencode_ = gencode.loc[gencode.gene_id == gene_id, bed_cols]
    snmf_meta = get_meta(snmf_)
    gencode_meta = get_meta(gencode_)

    snmf_unique_coords = get_difference_coords(snmf_meta, gencode_meta)
    gencode_unique_coords = get_difference_coords(gencode_meta, snmf_meta)

    

    if (len(snmf_unique_coords) >= 100) and (len(gencode_unique_coords) >= 100):

        
        counts = pd.read_csv(f'coverage/counts_filtered/{gene_id}.csv.gz', index_col=0)

        snmf_unique_coords = snmf_unique_coords.intersection(counts.columns)
        gencode_unique_coords = gencode_unique_coords.intersection(counts.columns)

        # snmf_unique_coords = get_counts_per_bp(snmf_unique_coords, counts)
        # gencode_unique_coords = get_counts_per_bp(gencode_unique_coords, counts)

        if (len(snmf_unique_coords) >= 100) and (len(gencode_unique_coords) >= 100):
            snmf_cbp = get_counts_per_bp(snmf_unique_coords, counts)
            gencode_cbp = get_counts_per_bp(gencode_unique_coords, counts)
            return snmf_cbp, gencode_cbp

def get_counts_per_bp_by_coords(counts, coord_array, min_bp):
    counts = np.array(counts).astype(int)
    counts = counts[coord_array].astype(int)
    len_counts = len(counts)
    if len_counts >= min_bp:
        cpm = np.sum(counts)/len_counts
    else:
        cpm = np.nan
    return cpm, len_counts

def get_residual_cpb(counts_fh, unique_coords_1, unique_coords_2, min_bp=100, first_group = 'snmf', second_group='gencode'):
    sample_list = []
    tissue_list = []
    cpm_coords_1_list = []
    cpm_coords_2_list = []
    len_coords_1_list = []
    len_coords_2_list = []
    first_line = True
    for line in counts_fh:
        line = np.array(line.decode().split(','))
        sample = line[0]
        line_counts = line[1:]
        if first_line:
            count_cols = pd.Index(line_counts)
            region_start = int(count_cols[0].split(':')[1])
            region_end = int(count_cols[-1].split(':')[1])
            region_length = region_end - region_start
            coord_array_1 = count_cols.isin(unique_coords_1)
            coord_array_2 = count_cols.isin(unique_coords_2)
            if (np.sum(coord_array_1) < min_bp) or  (np.sum(coord_array_2) < min_bp):
                return None
            first_line = False
        else:
            cpm_coords_1, len_coords_1 = get_counts_per_bp_by_coords(line_counts, coord_array_1, min_bp)
            cpm_coords_2, len_coords_2 = get_counts_per_bp_by_coords(line_counts, coord_array_2, min_bp)
            
            cpm_coords_1_list.append(cpm_coords_1)
            cpm_coords_2_list.append(cpm_coords_2)
            len_coords_1_list.append(len_coords_1)
            len_coords_2_list.append(len_coords_2)
            sample, tissue = sample.split('.')
            sample_list.append(sample)
            tissue_list.append(tissue)

    out_df = pd.DataFrame()
    out_df['samples'] = sample_list 
    out_df['tissue'] = tissue_list 
    out_df[first_group] = cpm_coords_1_list 
    out_df[second_group] = cpm_coords_2_list
    out_df['len_1'] = len_coords_1_list 
    out_df['len_2'] = len_coords_2_list
    out_df['gene_length'] = [region_length]*len(sample_list)
    # out_df['annotation'] = [first_group]*len(sample_list) + [second_group]*len(sample_list)

    return out_df
    

def cbp_wrapper(snmf_exons, gencode_exons, gene_id):
    bed_cols = ['chrom', 'start', 'end']
    snmf_meta_transcript = get_meta(snmf_exons.loc[snmf_exons.gene_id == gene_id, bed_cols])
    gencode_meta_transcript = get_meta(gencode_exons.loc[gencode_exons.gene_id == gene_id, bed_cols])
    snmf_unique_coords = get_difference_coords(snmf_meta_transcript, gencode_meta_transcript)
    gencode_unique_coords = get_difference_coords(gencode_meta_transcript, snmf_meta_transcript)
    if (len(snmf_unique_coords) >= 100) and (len(gencode_unique_coords) >= 100) and (len(snmf_unique_coords) <= 10000) and (len(gencode_unique_coords) <= 10000):
        counts_fh = gzip.open(f'/project2/mstephens/cfbuenabadn/gtex-stm/code/coverage/counts_filtered/{gene_id}.csv.gz')
        df = get_residual_cpb(counts_fh, snmf_unique_coords, gencode_unique_coords, min_bp=100)
        if df is None:
            return None
        df['gene_id'] = [gene_id]*len(df)
    else:
        df = None
    return df

if __name__ == '__main__':

    arguments = sys.argv
    gene_id = arguments[1]
    print(gene_id)

    gencode_exons_bed = '/project2/mstephens/cfbuenabadn/gtex-stm/code/gencode.v44.primary_assembly.exons.sorted.bed.gz'
    gencode_exons = pd.read_csv(gencode_exons_bed, sep='\t', names = ['chrom', 'start', 'end', 'gene_id', 
                                                           'transcript_id', 'strand', 'exon_id', 'transcript_support_level',
                                                           'basic', 'Ensembl_canonical', 'MANE_Select', 'appris'])
    
    snmf_exons_bed = '/project2/mstephens/cfbuenabadn/gtex-stm/code/snmf.exons.sorted.bed.gz'
    snmf_exons = pd.read_csv(snmf_exons_bed, sep='\t', names = ['chrom', 'start', 'end', 'gene_id', 
                                                     'transcript_id', 'strand', 'factors', 'exon_id'])

    
    annotations = ['primary_assembly', 'appris_principal', 'appris', 'basic', 'MANE_Select', 'Ensembl_canonical', 'support_level', 'basic']
    
    file_name = f'/project2/mstephens/cfbuenabadn/gtex-stm/code/ebpmf_models/filtered/unique_regions/K10/{annotation}/{gene_id}.unique_region_counts.tab.gz'

    with gzip.open(file_name, 'wb') as fh:
        first_line = True
        
        try:
            df = cbp_wrapper(snmf_exons, gencode_exons, gene_id)
            if df is None:
                print('no unique region of enough length')
            else:
                if first_line:
                    first_line_str = ('\t'.join(df.columns) + '\n').encode()
                    fh.write(first_line_str)
                    first_line = False
                    
                for idx, row in df.iterrows():
                    line = ('\t'.join(list(row.astype(str))) + '\n').encode()
                    fh.write(line)
        except:
            print('something went wrong')

    print('done!')
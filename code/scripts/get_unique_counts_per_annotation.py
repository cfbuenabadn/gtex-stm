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

def get_overlaps(exon_snmf, exon_appris, exon='first'):
    overlap = exon_snmf.intersect(exon_appris).to_dataframe()

    if len(overlap) == 0:
        return 0, 0
    
    print(overlap)
    print(len(overlap))
    overlap_bp = np.abs(int(overlap.iloc[0].start) - int(overlap.iloc[0].end))
    # print(overlap_bp)

    if exon == 'first':
        snmf_start = int(exon_snmf.to_dataframe().iloc[0].start)
        appris_start = int(exon_appris.to_dataframe().iloc[0].start)
        difference = np.abs(snmf_start - appris_start)
    elif exon == 'last':
        snmf_end = int(exon_snmf.to_dataframe().iloc[0].end)
        appris_end = int(exon_appris.to_dataframe().iloc[0].end)
        difference = np.abs(snmf_end - appris_end)
    else:
        raise Exception('wrong exon')

    return overlap_bp, difference
    
def correct_edges(snmf_meta, appris_meta):

    snmf_meta_out = snmf_meta.to_dataframe()
    appris_meta_out = appris_meta.to_dataframe()
    
    first_exon_snmf = BedTool.from_dataframe(pd.DataFrame(snmf_meta.to_dataframe().iloc[0]).T)
    first_exon_appris = BedTool.from_dataframe(pd.DataFrame(appris_meta.to_dataframe().iloc[0]).T)

    overlap_bp, difference = get_overlaps(first_exon_snmf, first_exon_appris, exon='first')

    if overlap_bp > 0:
        difference_ratio = difference/(difference + overlap_bp)

        # print(overlap_bp)
        # print(difference_ratio)

        if ((difference_ratio < 0.25) | (difference < 200)):
            max_start = int(np.max([int(first_exon_snmf.to_dataframe().iloc[0].start), int(first_exon_appris.to_dataframe().iloc[0].start)]))
            snmf_meta_out.loc[0, 'start'] = max_start
            appris_meta_out.loc[0, 'start'] = max_start
        elif (difference > 200):
            if int(first_exon_snmf.to_dataframe().iloc[0].start) > int(first_exon_appris.to_dataframe().iloc[0].start):
                snmf_meta_out.loc[0, 'start'] = int(snmf_meta_out.loc[0, 'start']) - 50
            else:
                appris_meta_out.loc[0, 'start'] = int(appris_meta_out.loc[0, 'start']) - 50
            

    last_exon_snmf = BedTool.from_dataframe(pd.DataFrame(snmf_meta.to_dataframe().iloc[len(snmf_meta.to_dataframe())-1]).T)
    last_exon_appris = BedTool.from_dataframe(pd.DataFrame(appris_meta.to_dataframe().iloc[len(appris_meta.to_dataframe())-1]).T)

    overlap_bp, difference = get_overlaps(last_exon_snmf, last_exon_appris, exon='last')

    if overlap_bp > 0:
        difference_ratio = difference/(difference + overlap_bp)

        if ((difference_ratio < 0.25) | (difference < 200)):
            min_end = int(np.max([int(last_exon_snmf.to_dataframe().iloc[0].end), int(last_exon_appris.to_dataframe().iloc[0].end)]))
            snmf_meta_out.loc[len(snmf_meta.to_dataframe())-1, 'end'] = min_end
            appris_meta_out.loc[len(appris_meta.to_dataframe())-1, 'end'] = min_end
        elif (difference > 200):
            if int(last_exon_snmf.to_dataframe().iloc[0].end) < int(last_exon_appris.to_dataframe().iloc[0].end):
                snmf_meta_out.loc[len(snmf_meta.to_dataframe())-1, 'end'] = int(snmf_meta_out.loc[len(snmf_meta.to_dataframe())-1, 'end']) + 50
            else:
                appris_meta_out.loc[len(appris_meta.to_dataframe())-1, 'end'] = int(appris_meta_out.loc[len(appris_meta.to_dataframe())-1, 'end']) + 50

    snmf_meta_out = BedTool.from_dataframe(snmf_meta_out)
    appris_meta_out = BedTool.from_dataframe(appris_meta_out)

    return snmf_meta_out, appris_meta_out


def get_difference_coords(gene1_bed, gene2_bed, min_diff = 50):
    sub_bed = gene1_bed.subtract(gene2_bed).to_dataframe()
    if len(sub_bed) > 0:
        print(sub_bed)
        sub_bed = sub_bed.loc[(sub_bed.end - sub_bed.start) > min_diff]
    sub_coords = get_transcript_coords(sub_bed)
    return sub_coords

def get_counts_per_bp(coords, counts):
    coords = counts.columns.intersection(coords)
    cpb = counts[coords].sum(axis=1)/len(coords)
    return cpb


def get_counts_per_bp_by_coords(counts, coord_array, min_bp):
    counts = np.array(counts).astype(int)
    counts = counts[coord_array].astype(int)
    len_counts = len(counts)
    if len_counts > 0:
        total_counts = np.sum(counts)
        cpm = total_counts/len_counts
    else:
        total_counts = 0
        cpm = np.nan
    return total_counts, cpm, len_counts

def get_residual_cpb(counts_fh, snmf_meta_transcript_per_run, gencode_metas_per_annotation, min_bp=100, first_group = 'snmf', second_group='gencode'):
    sample_list = []
    tissue_list = []
    counts_coords_1_list = []
    counts_coords_2_list = []
    cpm_coords_1_list = []
    cpm_coords_2_list = []
    len_coords_1_list = []
    len_coords_2_list = []
    annotation_list = []
    run_list = []
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
            #coord_array_1 = count_cols.isin(unique_coords_1)
            #coord_array_2 = count_cols.isin(unique_coords_2)
            # if (np.sum(coord_array_1) < min_bp) or (np.sum(coord_array_2) < min_bp):
            #     return None
            first_line = False
        else:
            sample, tissue = sample.split('.')
            for run in snmf_meta_transcript_per_run.keys():
                for annotation in gencode_metas_per_annotation.keys():
                    snmf_meta_transcript = snmf_meta_transcript_per_run[run]
                    gencode_meta_transcript = gencode_metas_per_annotation[annotation]
                    if (len(snmf_meta_transcript.to_dataframe()) > 0) and (len(gencode_meta_transcript.to_dataframe()) > 0):
                        snmf_meta_transcript, gencode_meta_transcript = correct_edges(snmf_meta_transcript, gencode_meta_transcript)
                    snmf_unique_coords = get_difference_coords(snmf_meta_transcript, gencode_meta_transcript)
                    gencode_unique_coords = get_difference_coords(gencode_meta_transcript, snmf_meta_transcript)
    
                    coord_array_1 = count_cols.isin(snmf_unique_coords)
                    coord_array_2 = count_cols.isin(gencode_unique_coords)
    
                    if (np.sum(coord_array_1) < min_bp) and (np.sum(coord_array_2) < min_bp):
                        continue
                    else:
                    
                        counts_coords_1, cpm_coords_1, len_coords_1 = get_counts_per_bp_by_coords(line_counts, coord_array_1, min_bp)
                        counts_coords_2, cpm_coords_2, len_coords_2 = get_counts_per_bp_by_coords(line_counts, coord_array_2, min_bp)
                        
                        counts_coords_1_list.append(counts_coords_1)
                        counts_coords_2_list.append(counts_coords_2)
                        cpm_coords_1_list.append(cpm_coords_1)
                        cpm_coords_2_list.append(cpm_coords_2)
                        len_coords_1_list.append(len_coords_1)
                        len_coords_2_list.append(len_coords_2)
                       
                        #print(sample)
                        #sample, tissue = sample.split('.')
                        sample_list.append(sample)
                        tissue_list.append(tissue)
                        annotation_list.append(annotation)
                        run_list.append(run)

    out_df = pd.DataFrame()
    out_df['samples'] = sample_list 
    out_df['tissue'] = tissue_list 
    out_df[first_group + '_cpb'] = cpm_coords_1_list 
    out_df[second_group + '_cpb'] = cpm_coords_2_list
    out_df[first_group + '_counts'] = counts_coords_1_list 
    out_df[second_group + '_counts'] = counts_coords_2_list
    out_df['len_1'] = len_coords_1_list 
    out_df['len_2'] = len_coords_2_list
    out_df['gene_length'] = [region_length]*len(sample_list)
    out_df['gencode_annotation'] = annotation_list
    out_df['snmf_annotation'] = run_list
    
    # out_df['annotation'] = [first_group]*len(sample_list) + [second_group]*len(sample_list)

    return out_df
    

def cbp_wrapper(snmf_exons, #snmf_exons_5, snmf_exons_3, 
                gencode_exons, gene_id):
    bed_cols = ['chrom', 'start', 'end']
    gencode_cols = ['#chrom', 'start', 'end']
    #annotations = ['primary_assembly', 'appris_principal', 'appris', 'basic', 'MANE_Select', 'Ensembl_canonical', 'support_level']

    snmf_meta_transcript = get_meta(snmf_exons.loc[snmf_exons.gene_id == gene_id, bed_cols])
    # snmf_meta_transcript_5 = get_meta(snmf_exons_5.loc[snmf_exons_5.gene_id == gene_id, bed_cols])
    # snmf_meta_transcript_3 = get_meta(snmf_exons_3.loc[snmf_exons_3.gene_id == gene_id, bed_cols])

    snmf_meta_transcript_per_run = {
        'snmf_10':snmf_meta_transcript#,
        # 'snmf_5':snmf_meta_transcript_5,
        # 'snmf_3':snmf_meta_transcript_3
    }
    
    gencode_transcript = gencode_exons.loc[gencode_exons.gene_id == gene_id]

    gencode_metas_per_annotation = {
        'primary_assembly':get_meta(gencode_transcript[gencode_cols]),
        'appris_principal':get_meta(gencode_transcript.loc[gencode_exons.appris.astype(str).apply(lambda x: 'appris_principal' in x), gencode_cols]),
        'appris':get_meta(gencode_transcript.loc[gencode_exons.appris.astype(str).apply(lambda x: 'appris' in x), gencode_cols]),
        'basic':get_meta(gencode_transcript.loc[gencode_exons.basic == 'basic', gencode_cols]),
        'MANE_Select':get_meta(gencode_transcript.loc[gencode_exons.MANE_Select == 'MANE_Select', gencode_cols]),
        'Ensembl_canonical':get_meta(gencode_transcript.loc[gencode_exons.Ensembl_canonical == 'Ensembl_canonical', gencode_cols]),
        'transcript_support_level':get_meta(gencode_transcript.loc[gencode_transcript.transcript_support_level > 0, gencode_cols])
    }
    
    
    #gencode_meta_transcript = get_meta(gencode_exons.loc[gencode_exons.gene_id == gene_id, bed_cols])
    #snmf_unique_coords = get_difference_coords(snmf_meta_transcript, gencode_meta_transcript)
    #gencode_unique_coords = get_difference_coords(gencode_meta_transcript, snmf_meta_transcript)
    #if (len(snmf_unique_coords) >= 50) and (len(gencode_unique_coords) >= 50) and (len(snmf_unique_coords) <= 100000) and (len(gencode_unique_coords) <= 100000):
    counts_fh = gzip.open(f'/project2/mstephens/cfbuenabadn/gtex-stm/code/coverage/counts_filtered/{gene_id}.csv.gz')
    df = get_residual_cpb(counts_fh, snmf_meta_transcript_per_run, gencode_metas_per_annotation, min_bp=50)
    if df is None:
        return None
    df['gene_id'] = [gene_id]*len(df)
    # else:
    #     df = None
    return df

if __name__ == '__main__':

    arguments = sys.argv
    gene_id = arguments[1]
    print(gene_id)

    gencode_exons_bed = '/project2/mstephens/cfbuenabadn/gtex-stm/code/gencode.v44.primary_assembly.exons.sorted.bed.gz'
    gencode_exons = pd.read_csv(gencode_exons_bed, sep='\t')#, names = ['chrom', 'start', 'end', 'gene_id', 
                                                           # 'transcript_id', 'strand', 'exon_id', 'transcript_support_level',
                                                           # 'basic', 'Ensembl_canonical', 'MANE_Select', 'appris'])
    
    snmf_exons_bed = '/project2/mstephens/cfbuenabadn/gtex-stm/code/ebpmf_models/filtered/snmf_10/tables/snmf.exons.sorted.bed.gz'
    snmf_exons = pd.read_csv(snmf_exons_bed, sep='\t', names = ['chrom', 'start', 'end', 'gene_id', 
                                                     'transcript_id', 'strand', 'factors', 'exon_id'])

    # snmf_exons_bed5 = '/project2/mstephens/cfbuenabadn/gtex-stm/code/snmf_5.exons.sorted.bed.gz'
    # snmf_exons_5 = pd.read_csv(snmf_exons_bed5, sep='\t', names = ['chrom', 'start', 'end', 'gene_id', 
    #                                                  'transcript_id', 'strand', 'factors', 'exon_id'])

    # snmf_exons_bed3 = '/project2/mstephens/cfbuenabadn/gtex-stm/code/snmf_3.exons.sorted.bed.gz'
    # snmf_exons_3 = pd.read_csv(snmf_exons_bed3, sep='\t', names = ['chrom', 'start', 'end', 'gene_id', 
    #                                                  'transcript_id', 'strand', 'factors', 'exon_id'])

    
    
    
    file_name = f'/project2/mstephens/cfbuenabadn/gtex-stm/code/ebpmf_models/filtered/snmf_10/unique_regions/{gene_id}.unique_region_counts.tab.gz'

    with gzip.open(file_name, 'wb') as fh:
        first_line = True
        
        # try:
        df = cbp_wrapper(snmf_exons, #snmf_exons_5, snmf_exons_3, 
                         gencode_exons, gene_id)
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
        # except:
        #     print('something went wrong')

    print('done!')

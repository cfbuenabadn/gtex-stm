import numpy as np
import pandas as pd
import gzip
from pybedtools import BedTool
import os
from scipy.stats import zscore


def aggregate_dataframe(df):
    grouped_df = df.groupby(['chrom', 'start', 'end'])[['chrom', 'start', 'end', 'chrom_exon', 'start_exon', 'end_exon']].apply(lambda group: pd.Series({
        'chrom_exon_1': group.loc[group['start_exon'].idxmin(), 'chrom_exon'],
        'start_exon_1': group.loc[group['start_exon'].idxmin(), 'start_exon'],
        'end_exon_1': group.loc[group['start_exon'].idxmin(), 'end_exon'],
        'chrom_exon_2': group.loc[group['start_exon'].idxmax(), 'chrom_exon'],
        'start_exon_2': group.loc[group['start_exon'].idxmax(), 'start_exon'],
        'end_exon_2': group.loc[group['start_exon'].idxmax(), 'end_exon'],
    })).reset_index()

    grouped_df.start += 2
    grouped_df.end -= 2
    
    return grouped_df

def get_gencode_bed(gencode, gene):
    gencode_sub = gencode.loc[(gencode.gene == gene)].dropna()
    if gencode_sub.appris.apply(lambda x: 'appris_principal_' in x).any():
        gencode_transcript = gencode_sub.loc[gencode_sub.appris.apply(lambda x: 'appris_principal_' in x)]
    else:
        gencode_transcript = gencode_sub.loc[gencode_sub.mane == 'MANE_Select']

    gencode_bed = BedTool.from_dataframe(gencode_transcript).merge()

    gencode_full_gene = pd.DataFrame()
    gencode_full_gene['chrom'] = [gencode_transcript.iloc[0, 0]]
    gencode_full_gene['start'] = [gencode_transcript.start.min()]
    gencode_full_gene['end'] = [gencode_transcript.end.max()]

    gencode_full_gene = BedTool.from_dataframe(gencode_full_gene)
    
    return gencode_bed, gencode_full_gene

def get_IR_bed(torino_bed, gencode_bed, gencode_full_gene):
    subtracted = torino_bed.intersect(gencode_full_gene, f=1).subtract(gencode_bed).to_dataframe()
    subtracted.start -= 1
    subtracted.end += 1
    df_ = BedTool.from_dataframe(subtracted).intersect(gencode_bed, c=True).to_dataframe(
        names = ['chrom', 'start', 'end', 'gene', 'isoform', 'strand', 'factors', 'exon', 'overlaps']
    )
    
    overlaps_2 = df_.loc[df_.overlaps == 2]
    df_ = BedTool.from_dataframe(overlaps_2).intersect(gencode_bed, wo=True).to_dataframe(
    names =['chrom', 'start', 'end', 'gene', 'isoform', 'strand', 'factors', 'exon', 'overlaps', 'chrom_exon', 'start_exon', 'end_exon', 'overlap_bp'])

    agg_dataframe = aggregate_dataframe(df_)
    return agg_dataframe


def get_counts_for_intervals_grouped(df, coords, counts):
    result = []
    
    # Iterate through each interval in the dataframe
    for _, row in df.iterrows():
        start, end = row['start'], row['end']
        start_exon_1, end_exon_1 = row['start_exon_1'], row['end_exon_1']
        start_exon_2, end_exon_2 = row['start_exon_2'], row['end_exon_2']
        
        # Get indices of coordinates that fall within the interval
        in_interval_idx = [i for i, coord in enumerate(coords) if start <= coord <= end]
        # Get the corresponding counts for those coordinates
        in_interval_counts = [counts[i] for i in in_interval_idx]
        
        # Determine the flanking intervals using start_exon_1/end_exon_1 and start_exon_2/end_exon_2
        flanking_before_idx = [i for i, coord in enumerate(coords) if start_exon_1 <= coord <= end_exon_1]
        flanking_after_idx = [i for i, coord in enumerate(coords) if start_exon_2 <= coord <= end_exon_2]
        
        flanking_before = [counts[i] for i in flanking_before_idx]
        flanking_after = [counts[i] for i in flanking_after_idx]
        
        result.append({
            'interval': (start, end),
            'in_interval': in_interval_counts,
            'flanking_before': flanking_before,
            'flanking_after': flanking_after
        })
    
    return result


def get_plotting_averages2(isoform_name, aggregated, metadata, coords_int, fh):
    
    intron_idx = get_intron_idx(coords_int, aggregated, isoform_name)

    samples_list = []
    tissue_list = []
    AD_list = []
    braak_list = []
    condition_list = []
    tissue_list = []
    for line in fh:
        sample_and_counts = line.decode().rstrip().split(',')
        sample = sample_and_counts[0]
        if sample not in metadata.index:
            continue
        
        counts_sample = np.array(sample_and_counts[1:]).astype(int)
        sample_meta = metadata.loc[sample]
        braaksc = metadata.loc[sample].loc['braaksc']
        ad = metadata.loc[sample].loc['AD']
        condition = metadata.loc[sample].loc['condition']
        tissue = metadata.loc[sample].loc['tissue_id']

        samples_list.append(sample)
        braak_list.append(braaksc)
        AD_list.append(ad)
        condition_list.append(condition)
        tissue_list.append(tissue)

        for intron, val in intron_idx.items():
            in_interval = val['intron']['idx']
            len_intron = val['intron']['len']
            counts_in_interval = np.sum(counts_sample[in_interval])
            coverage_intron = counts_in_interval/len_intron

            in_flk1 = val['flk1']['idx']
            in_flk2 = val['flk2']['idx']
            len_flk1 = val['flk1']['len']
            len_flk2 = val['flk2']['len']

            counts_in_exons = np.sum(counts_sample[in_flk1]) + np.sum(counts_sample[in_flk2])
            
            coverage_exons = counts_in_exons/(len_flk1 + len_flk2)
            
            if (coverage_exons == 0):
                if (coverage_intron > 1e-4):
                    PSI = 1
                else:
                    PSI = 0
            else:
                PSI = np.min([1, coverage_intron/coverage_exons])

            intron_idx[intron]['coverage'].append(PSI)

    meta_samples = pd.DataFrame()
    meta_samples['sample'] = samples_list
    meta_samples['braaksc'] = braak_list
    meta_samples['AD'] = AD_list
    meta_samples['condition'] = condition_list
    meta_samples['tissue'] = tissue_list
    return intron_idx, meta_samples

    
def get_intron_idx(coords_int, aggregated, isoform_name):
    introns_dict = {}
   
        
    for idx, row in aggregated.iterrows():
        interval = (row.chrom, row.start, row.end, isoform_name.split('.')[0])
        idx_intron = (np.array(coords_int) >= row.start) & (np.array(coords_int) <= row.end)
        idx_flk1 = (np.array(coords_int) >= row.start_exon_1) & (np.array(coords_int) <= row.end_exon_1)
        idx_flk2 = (np.array(coords_int) >= row.start_exon_2) & (np.array(coords_int) <= row.end_exon_2)
        intron = isoform_name + f'.{str(idx+1)}'
        intron_interval = {'idx':idx_intron, 'len':np.sum(idx_intron)}
        flk1_interval = {'idx':idx_flk1, 'len':np.sum(idx_flk1)}
        flk2_interval = {'idx':idx_flk2, 'len':np.sum(idx_flk2)}
        introns_dict.update({intron:{'interval':interval, 'intron':intron_interval, 'flk1':flk1_interval, 'flk2':flk2_interval, 'coverage':[]}})
    return introns_dict




def prepare_intron_coverage(coverage, condition_list, intron):
    # coverage = np.array(intron_idx[intron]['coverage'])
    zcoverage = zscore(coverage)

    coverage_list = []
    for lista in condition_list:
        coverage_list.append(coverage[lista].mean())
        coverage_list.append(zcoverage[lista].mean())

    out_list = [str(x) for x in coverage_list]
    return out_list
    
def prepare_gene_df(intron_idx, meta_samples):
    braak_list = [list(meta_samples.braaksc == 0), 
                      list(meta_samples.braaksc == 1), 
                      list(meta_samples.braaksc == 2), 
                      list(meta_samples.braaksc == 3),
                      list(meta_samples.braaksc == 4),
                      list(meta_samples.braaksc == 5),
                      list(meta_samples.braaksc == 6),
                      list(meta_samples.braaksc <= 1),
                      list(meta_samples.braaksc >= 5)]

    meta_samples_ac = meta_samples.loc[meta_samples.tissue=='AC']
    braak_ac_list = [list(meta_samples_ac.braaksc == 0), 
                      list(meta_samples_ac.braaksc == 1), 
                      list(meta_samples_ac.braaksc == 2), 
                      list(meta_samples_ac.braaksc == 3),
                      list(meta_samples_ac.braaksc == 4),
                      list(meta_samples_ac.braaksc == 5),
                      list(meta_samples_ac.braaksc == 6),
                      list(meta_samples_ac.braaksc <= 1),
                      list(meta_samples_ac.braaksc >= 5)]

    meta_samples_dlpc = meta_samples.loc[meta_samples.tissue=='DLPC']
    braak_dlpc_list = [list(meta_samples_dlpc.braaksc == 0), 
                      list(meta_samples_dlpc.braaksc == 1), 
                      list(meta_samples_dlpc.braaksc == 2), 
                      list(meta_samples_dlpc.braaksc == 3),
                      list(meta_samples_dlpc.braaksc == 4),
                      list(meta_samples_dlpc.braaksc == 5),
                      list(meta_samples_dlpc.braaksc == 6),
                      list(meta_samples_dlpc.braaksc <= 1),
                      list(meta_samples_dlpc.braaksc >= 5)]

    meta_samples_pcc = meta_samples.loc[meta_samples.tissue=='PCC']
    braak_pcc_list = [list(meta_samples_pcc.braaksc == 0), 
                      list(meta_samples_pcc.braaksc == 1), 
                      list(meta_samples_pcc.braaksc == 2), 
                      list(meta_samples_pcc.braaksc == 3),
                      list(meta_samples_pcc.braaksc == 4),
                      list(meta_samples_pcc.braaksc == 5),
                      list(meta_samples_pcc.braaksc == 6),
                      list(meta_samples_pcc.braaksc <= 1),
                      list(meta_samples_pcc.braaksc >= 5)]
    
    ad_list = [list(meta_samples.condition == 'no condition'),
               list(meta_samples.condition == 'AD')]
    ad_ac_list = [list(meta_samples_ac.condition == 'no condition'),
                  list(meta_samples_ac.condition == 'AD')]
    ad_dlpc_list = [list(meta_samples_dlpc.condition == 'no condition'),
                    list(meta_samples_dlpc.condition == 'AD')]
    ad_pcc_list = [list(meta_samples_pcc.condition == 'no condition'),
                   list(meta_samples_pcc.condition == 'AD')]

    
    list_of_rows = []
    for intron, intron_data in intron_idx.items():

        coverage = np.array(intron_idx[intron]['coverage'])
        coverage_ac = np.array(intron_idx[intron]['coverage'])[list(meta_samples.tissue=='AC')]
        coverage_dlpc = np.array(intron_idx[intron]['coverage'])[list(meta_samples.tissue=='DLPC')]
        coverage_pcc = np.array(intron_idx[intron]['coverage'])[list(meta_samples.tissue=='PCC')]

        braak_coverage_list = prepare_intron_coverage(coverage, braak_list, intron)
        braak_ac_coverage_list = prepare_intron_coverage(coverage_ac, braak_ac_list, intron)
        braak_dlpc_coverage_list = prepare_intron_coverage(coverage_dlpc, braak_dlpc_list, intron)
        braak_pcc_coverage_list = prepare_intron_coverage(coverage_pcc, braak_pcc_list, intron)

        ad_coverage_list =  prepare_intron_coverage(coverage, ad_list, intron)
        ad_ac_coverage_list = prepare_intron_coverage(coverage_ac, ad_ac_list, intron)
        ad_dlpc_coverage_list = prepare_intron_coverage(coverage_dlpc, ad_dlpc_list, intron)
        ad_pcc_coverage_list = prepare_intron_coverage(coverage_pcc, ad_pcc_list, intron)
        
        coverage_list = braak_coverage_list + braak_ac_coverage_list + braak_dlpc_coverage_list + braak_pcc_coverage_list
        coverage_list += ad_coverage_list + ad_ac_coverage_list + ad_dlpc_coverage_list + ad_pcc_coverage_list
        #prepare_intron_coverage(intron_idx, condition_list, intron)
        
        bed_file = list(intron_data['interval'])
        bed_file += ['.'.join(intron.split('.')[:2])]
        bed_file += coverage_list
        bed_file = [str(x) for x in bed_file]
        bed_row = '\t'.join(bed_file) + '\n'
        bed_row = bed_row.encode()
        list_of_rows.append(bed_row)

    return list_of_rows
        


isoforms = pd.read_csv('ebpmf_models/gao_models/tables/snmf_5.merged_isoforms.exons.sorted.bed.gz', sep='\t',
                      names = ['chrom', 'start', 'end', 'gene_id', 'transcript_id', 'strand', 'factor', 'exon'])

gencode = pd.read_csv('Annotations/gencode.v44.primary_assembly.exons.sorted.bed.gz', sep='\t', 
                     names = ['chrom', 'start', 'end', 'gene', 'transcript', 'strand', 'exon', 'score', 'gencode_annot', 
                              'ensembl', 'mane', 'appris', 'annotation'])

ad_boxplot = pd.read_csv('intron_retention.second_AD_dataframe.tab.gz', sep='\t')

EL = pd.read_csv('../code/ebpmf_models/gao_models/tables/snmf5_EL.tab.gz', sep='\t', index_col=3)


samples = EL.columns[5:]

metadata = pd.read_csv('../code/gao_data/metadata/merged_metadata.tab.gz', sep='\t', index_col=0)

X = EL[EL.columns[5:]]

metadata_test = metadata.loc[EL.columns.intersection(metadata.index)]


column_names = ['chrom', 'start', 'end', 'gene', 'isoform']
braak_columns = []
for i in range(7):
    braaks_col = f'braaks_{str(i)}_loadings'
    braak_columns.append(braaks_col)
    braaks_col = f'braaks_{str(i)}_zloadings'
    braak_columns.append(braaks_col)

braak_columns.extend(['braaks_low_loadings',
                       'braaks_low_zloadings',
                       'braaks_high_loadings',
                       'braaks_high_zloadings'])

braak_ac_columns = [x + '.AC' for x in braak_columns]
braak_dlpc_columns = [x + '.DLPC' for x in braak_columns]
braak_pcc_columns = [x + '.PCC' for x in braak_columns]

ad_columns = ['ctrl_loadings',
              'ctrl_zloadings',
              'ad_loadings',
              'ad_zloadings']

ad_ac_columns = [x + '.AC' for x in ad_columns]
ad_dlpc_columns = [x + '.DLPC' for x in ad_columns]
ad_pcc_columns = [x + '.PCC' for x in ad_columns]

column_names += braak_columns + braak_ac_columns + braak_dlpc_columns + braak_pcc_columns + ad_columns + ad_ac_columns + ad_dlpc_columns + ad_pcc_columns
    

import logging

logging.basicConfig(filename='logs/print_second_ad_introns_bed.log', level=logging.INFO, 
                    format='%(asctime)s - %(message)s')



with gzip.open('ebpmf_models/gao_models/tables/second_AD_IR_plots_1e3.bed.gz', 'wb') as fh2:

    column_row = ('\t'.join(column_names) + '\n').encode()
    fh2.write(column_row)

    for idx, df in ad_boxplot.loc[ad_boxplot.FDR <= 1e-3].groupby(['gene', 'isoform']):
        try:
            gene, isoform = idx
            gencode_bed, gencode_full_gene = get_gencode_bed(gencode, gene)
            torino_bed = BedTool.from_dataframe(isoforms.loc[isoforms.transcript_id == isoform])
            IR_bed = get_IR_bed(torino_bed, gencode_bed, gencode_full_gene)
    
            fh = gzip.open(f'gao_data/counts/{gene}.csv.gz')
            coords = fh.readline()
            coords_int = [int(x.split(':')[1]) for x in coords.decode().rstrip().split(',')[1:]]
            intron_idx, meta_samples = get_plotting_averages2(isoform, IR_bed, metadata_test, coords_int, fh)
            list_of_rows = prepare_gene_df(intron_idx, meta_samples)
            fh2.writelines(list_of_rows)
            fh.close()
            logging.info(f"processed {gene}")
        except:
            continue

logging.info(f"Finished 1e-3 genes")

with gzip.open('ebpmf_models/gao_models/tables/second_AD_IR_plots_1e2.bed.gz', 'wb') as fh2:

    column_row = ('\t'.join(column_names) + '\n').encode()
    fh2.write(column_row)

    for idx, df in ad_boxplot.loc[(ad_boxplot.FDR > 1e-3) & (ad_boxplot.FDR <= 1e-2)].groupby(['gene', 'isoform']):
        try:
            gene, isoform = idx
            gencode_bed, gencode_full_gene = get_gencode_bed(gencode, gene)
            torino_bed = BedTool.from_dataframe(isoforms.loc[isoforms.transcript_id == isoform])
            IR_bed = get_IR_bed(torino_bed, gencode_bed, gencode_full_gene)
    
            fh = gzip.open(f'gao_data/counts/{gene}.csv.gz')
            coords = fh.readline()
            coords_int = [int(x.split(':')[1]) for x in coords.decode().rstrip().split(',')[1:]]
            intron_idx, meta_samples = get_plotting_averages2(isoform, IR_bed, metadata_test, coords_int, fh)
            list_of_rows = prepare_gene_df(intron_idx, meta_samples)
            fh2.writelines(list_of_rows)
            fh.close()
            logging.info(f"processed {gene}")
        except:
            continue

logging.info(f"Finished 1e-2 genes")


with gzip.open('ebpmf_models/gao_models/tables/second_AD_IR_plots_5e2.bed.gz', 'wb') as fh2:

    column_row = ('\t'.join(column_names) + '\n').encode()
    fh2.write(column_row)

    for idx, df in ad_boxplot.loc[(ad_boxplot.FDR > 1e-2)].groupby(['gene', 'isoform']):
        try:
            gene, isoform = idx
            gencode_bed, gencode_full_gene = get_gencode_bed(gencode, gene)
            torino_bed = BedTool.from_dataframe(isoforms.loc[isoforms.transcript_id == isoform])
            IR_bed = get_IR_bed(torino_bed, gencode_bed, gencode_full_gene)
    
            fh = gzip.open(f'gao_data/counts/{gene}.csv.gz')
            coords = fh.readline()
            coords_int = [int(x.split(':')[1]) for x in coords.decode().rstrip().split(',')[1:]]
            intron_idx, meta_samples = get_plotting_averages2(isoform, IR_bed, metadata_test, coords_int, fh)
            list_of_rows = prepare_gene_df(intron_idx, meta_samples)
            fh2.writelines(list_of_rows)
            fh.close()
            logging.info(f"processed {gene}")
        except:
            continue

logging.info(f"Finished 5e-2 genes")

    




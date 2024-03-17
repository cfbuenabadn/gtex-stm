import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import os
import subprocess as sp
from tqdm import tqdm 
import gtfparse

import seaborn as sns

import tabix
from pybedtools import BedTool

def read_EF(gene, K):
    EF_smooth = pd.read_csv(f'ebpmf_models/tables/{gene}_K{str(K)}.tab.gz', sep='\t', index_col=0)
#     EF_smooth = pd.read_csv(f'ebpmf_models/tables/{gene}/ebpmf_K{str(K)}/EF_smooth.tab.gz', sep='\t', index_col=0)
#     col_names = [f'factor{str(i+1)}' for i in range(K)]
#     EF.columns = col_names
#     EF_smooth.columns = col_names

#     EF = EF.iloc[1:(len(EF)-1)]
#     EF_smooth = EF_smooth.iloc[1:(len(EF_smooth)-1)]
#     EF_smooth.index = EF.index
    
    step = int(EF_smooth.index[1].split(':')[1]) - int(EF_smooth.index[0].split(':')[1])
    
    if step > 1:
        EF_smooth = extend_EF(EF_smooth)
    
    step = int(EF_smooth.index[1].split(':')[1]) - int(EF_smooth.index[0].split(':')[1])
    
    if step > 1:
        EF_smooth = extend_EF(EF_smooth)
    
    return EF_smooth

# def read_EF(gene, K):
#     EF = pd.read_csv(f'ebpmf_models/tables/{gene}/ebpmf_K{str(K)}/EF.tab.gz', sep='\t', index_col=0)
#     EF_smooth = pd.read_csv(f'ebpmf_models/tables/{gene}/ebpmf_K{str(K)}/EF_smooth.tab.gz', sep='\t', index_col=0)
#     col_names = [f'factor{str(i+1)}' for i in range(K)]
#     EF.columns = col_names
#     EF_smooth.columns = col_names

#     EF = EF.iloc[1:(len(EF)-1)]
#     EF_smooth = EF_smooth.iloc[1:(len(EF_smooth)-1)]
#     EF_smooth.index = EF.index
    
#     step = int(EF_smooth.index[1].split(':')[1]) - int(EF_smooth.index[0].split(':')[1])
    
#     if step > 1:
#         EF_smooth = extend_EF(EF_smooth)
    
#     step = int(EF_smooth.index[1].split(':')[1]) - int(EF_smooth.index[0].split(':')[1])
    
#     if step > 1:
#         EF_smooth = extend_EF(EF_smooth)
    
#     return EF_smooth

# def extend_EF(EF):
#     print('extending')
#     step = int(EF.index[1].split(':')[1]) - int(EF.index[0].split(':')[1])
    
#     extended_idx = []
#     factor1 = []
#     factor2 = []
    
#     for idx, row in EF.iterrows():
#         chrom = idx.split(':')[0]
#         start = int(idx.split(':')[1])
#         extended_idx_ = [f'{chrom}:{str(i)}' for i in range(start, start+step)]
#         factor1_ = [row.factor1]*step
#         factor2_ = [row.factor2]*step
        
#         extended_idx.extend(extended_idx_)
#         factor1.extend(factor1_)
#         factor2.extend(factor2_)
        
#     extended_EF = pd.DataFrame()
#     extended_EF['factor1'] = factor1
#     extended_EF['factor2'] = factor2
#     extended_EF.index = extended_idx
    
#     return extended_EF


def extend_EF(EF):
    step = int(EF.index[1].split(':')[1]) - int(EF.index[0].split(':')[1])
    
    extended_idx = []
    factor1 = []
    factor2 = []
    factor3 = []
    factor4 = []
    factor5 = []
    
    for idx, row in EF.iterrows():
        chrom = idx.split(':')[0]
        start = int(idx.split(':')[1])
        extended_idx_ = [f'{chrom}:{str(i)}' for i in range(start, start+step)]
        factor1_ = [row.factor1]*step
        factor2_ = [row.factor2]*step
        factor3_ = [row.factor3]*step
        factor4_ = [row.factor4]*step
        factor5_ = [row.factor5]*step
        
        extended_idx.extend(extended_idx_)
        factor1.extend(factor1_)
        factor2.extend(factor2_)
        factor3.extend(factor3_)
        factor4.extend(factor4_)
        factor5.extend(factor5_)
        
    extended_EF = pd.DataFrame()
    extended_EF['factor1'] = factor1
    extended_EF['factor2'] = factor2
    extended_EF['factor3'] = factor3
    extended_EF['factor4'] = factor4
    extended_EF['factor5'] = factor5
    extended_EF.index = extended_idx
    
    return extended_EF


def binarize_EF(EF):
    
    EF_smooth = pd.DataFrame()
    for factor in EF.columns:
        f_smooth = smooth_factor(EF[factor])
        EF_smooth[factor] = f_smooth
        
    EF_smooth.index = EF.index
    
    
    return EF_smooth


########################
# def get_factor_shift(gene, K):
#     EF = read_EF(gene, K)
#     EF_smooth = smooth_EF(EF)
#     EF_shift = (EF_smooth - EF_smooth.shift(-1)).fillna(0)
#     return EF, EF_smooth, EF_shift
# 
###############################

def get_factor_shift(EF_smooth):
    EF_shift = (EF_smooth - EF_smooth.shift(-1)).fillna(0)
    return EF_shift


def make_factor_junctions_bed(EF_shift, factor='factor1'):
    
    EF_shift = EF_shift.loc[EF_shift[factor] != 0]
    
    open_ = False
    bed_factor = pd.DataFrame()
    chrom_list = []
    start_list = []
    end_list = []

    for idx, row in EF_shift.iterrows():
        if row[factor] == 1:
            chrom_, start_ = idx.split(':')
            start_ = int(start_)
            open_ = True
            chrom_list.append(chrom_)
            start_list.append(start_)
        elif (row[factor] == -1) and open_:
            chrom_, end_ = idx.split(':')
            end_ = int(end_)
            end_list.append(end_)

    factor_junctions_bed = pd.DataFrame([chrom_list, start_list, end_list]).T.dropna()
    factor_junctions_bed = BedTool.from_dataframe(factor_junctions_bed)
    
    return factor_junctions_bed

def run_tabix(factor, window_len = 50):
        
    chrom_ = factor.chrom
    start_ = int(factor.start) - window_len
    end_ = int(factor.start) + window_len

    tb = tabix.open('/project2/mstephens/cfbuenabadn/gtex-stm/code/junctions.tab.gz')
    juncs = tb.query(chrom_, start_, end_)
    
    return juncs

def get_overlapping_junctions_bed(factor):
    
    juncs = run_tabix(factor)
    
    df = pd.DataFrame()

    for idx, record in enumerate(juncs):
        df[f'record{str(idx)}'] = record[:-1]

    df = df.T
    overlapping_junctions_bed = BedTool.from_dataframe(df)
    
    return overlapping_junctions_bed

def make_junction_directory(idx, factor_junction, matching_junction):
    distance_start = np.abs(int(matching_junction.start) - int(factor_junction.start))
    distance_end = np.abs(int(matching_junction.end) - int(factor_junction.end))
    total_distance = distance_start + distance_end

    median_counts = np.quantile([int(x) for x in matching_junction.fields[4:]], 0.95)

    junction_dict_idx = {idx:{
        'start':int(matching_junction.start),
        'end':int(matching_junction.end),
        'distance':total_distance,
        'counts':median_counts
    }}
    
    return junction_dict_idx

def get_best_junction(factor_junction, intersection):
    
    junction_dict = {}
    
    for idx, matching_junction in enumerate(intersection):
        junction_dict_idx = make_junction_directory(idx, factor_junction, matching_junction)
        junction_dict.update(junction_dict_idx)
    
    current_junc = None
    for idx, junc in junction_dict.items():
        if (junc['distance'] < 100) and (junc['counts'] > 20):
            if current_junc:
                if (junc['distance'] < current_junc['distance']):
                    current_junc = junc
                else:
                    continue
            else:
                current_junc = junc
                
    return current_junc
        

def get_factor_junctions(factor_junctions_bed):
    
    junctions_match_list = []
    
    for factor_junction in factor_junctions_bed:
        
        overlapping_junctions_bed = get_overlapping_junctions_bed(factor_junction)
        intersection = overlapping_junctions_bed.intersect(factor_junctions_bed, f=0.9, r=True, wa=True)
        
        if len(intersection) == 0:
            junctions_match_list.append(None)
        else:
            best_junc = get_best_junction(factor_junction, intersection)
            junctions_match_list.append(best_junc)
            
    return junctions_match_list


def get_factor_bed(EF_shift, factor):
    open_ = False
    bed_factor = pd.DataFrame()
    chrom_list = []
    start_list = []
    end_list = []

    for idx, row in EF_shift.iterrows():
        if row[factor] == 1:
            chrom_, start_ = idx.split(':')
            start_ = int(start_)
            open_ = True
            chrom_list.append(chrom_)
            start_list.append(start_)
        elif (row[factor] == -1) and open_:
            chrom_, end_ = idx.split(':')
            end_ = int(end_)
            end_list.append(end_)

    bed_factors_df = pd.DataFrame([chrom_list, start_list, end_list]).T.dropna()
    bed_factors = BedTool.from_dataframe(bed_factors_df)
    
    return bed_factors_df, bed_factors

def correct_factor(EF, EF_smooth, EF_shift, factor, factor_junctions, bed_factors_df, cutoff_strict = 0.05, exon_quant = 0.99):
    
    EF_norm = EF/EF[factor].quantile(exon_quant)
    
    
    corrected_factor = []
    
    idx_covered = []
        
    for current_junc_pos, current_junc in enumerate(factor_junctions):
        
        factor_gap = bed_factors_df.iloc[current_junc_pos]
        chrom_gap = factor_gap[0]
        start_gap = int(factor_gap[1])
        end_gap = int(factor_gap[2])
        
        
        if current_junc:
            start_junc = current_junc['start']
            end_junc = current_junc['end']
            junc_length = end_junc - start_junc
            gap_coverage = [0]*junc_length
            
            junc_idx = [f'{chrom_gap}:{str(x)}' for x in range(start_junc, end_junc)]
            
        
        
        else:
            start_junc = start_gap
            end_junc = end_gap
            
            gap_length = end_gap - start_gap
            
            start_coverage_gap = start_gap - gap_length
            
            start_coverage = int(EF_shift.iloc[current_junc_pos].name.split(':')[1])
            
#             print(start_coverage)
#             print(start_gap)
            start_coverage = np.max([start_coverage_gap, start_coverage])
            
            exon_idx = [f'{chrom_gap}:{str(x)}' for x in range(start_coverage, start_gap)]
            gap_idx = [f'{chrom_gap}:{str(x)}' for x in range(start_gap, end_gap)]
            
            
            exon_median = EF_norm.loc[exon_idx, factor].median()
            intron_median = EF_norm.loc[gap_idx, factor].median()
            
            junc_idx = gap_idx
            
            if (intron_median > (exon_median*0.1)) & (intron_median > cutoff_strict):
                gap_coverage = [1]*len(gap_idx)
#                 print(list(EF_smooth.loc[gap_idx, factor]))
            else:
#                 print(gap_idx)
#                 print(type(gap_idx))
                gap_coverage = list(EF_smooth.loc[gap_idx, factor])
                
        if current_junc_pos == 0:
            _, first_pos = EF_smooth.index[0].split(':')
            first_pos = int(first_pos)
            idx_5p = [f'{chrom_gap}:{str(x)}' for x in range(first_pos, start_junc)]
            EF_5p = list(EF_smooth.loc[idx_5p, factor])
            for i, bp in reversed(list(enumerate(EF_5p))):
                if bp == 0:
                    EF_5p[i] = 1
                else:
                    break
            corrected_factor.extend(EF_5p)
            last_exon_5p = end_junc
            idx_covered.extend(idx_5p)
            

        else:
            len_exon_5p = start_junc - last_exon_5p
            exon_5p = [1]*len_exon_5p
            corrected_factor.extend(exon_5p)
            
            idx_5p = [f'{chrom_gap}:{str(x)}' for x in range(last_exon_5p, start_junc)]
            last_exon_5p = end_junc
            
            
            
            idx_covered.extend(idx_5p)

            
        idx_covered.extend(junc_idx)
            
        corrected_factor.extend(gap_coverage)

        
        
        if current_junc_pos == (len(factor_junctions)-1):
            _, last_pos = EF_smooth.index[-1].split(':')
            last_pos = int(last_pos)+1
            idx_3p = [f'{chrom_gap}:{str(x)}' for x in range(end_junc, last_pos)]
            EF_3p = list(EF_smooth.loc[idx_3p, factor])
            
            for i, bp in enumerate(EF_3p):
                if bp == 0:
                    EF_3p[i] = 1
                else:
                    break
            
            
            corrected_factor.extend(EF_3p)
            idx_covered.extend(idx_3p)

    corrected_factor = np.array(corrected_factor)
            
    return corrected_factor, idx_covered
            

def smooth_factor(factor, cutoff_=0.25, smooth_fraction = 0.25, step_fraction = 10, exon_quant = 0.99, cutoff_strict = 0.1):
    
    cutoff = cutoff_
    
    factor = factor/(np.quantile(factor, exon_quant))
    step =int(len(factor)/step_fraction)
        
    factor_smooth = []
    last_pos = None
    for idx, f in enumerate(factor):

        if f <= cutoff_strict:
            factor_smooth.append(0)
            last_pos = None
            continue
            
        if (idx == 0):
            if (f > cutoff):
                f_smooth = 1
                last_pos = 0
            else:
                f_smooth = 0
                last_pos = None
                
            
        else:
            
            if (factor_smooth[idx-1] == 0):
                if (f < cutoff):
                    f_smooth = 0
                    last_pos = None
                else:
                    f_smooth = 1
                    last_pos = idx
            else:
                last_pos = np.max([last_pos, (idx-100)])
                f_exon_peak = np.quantile(factor[last_pos:idx], exon_quant)
                
                if (f < (f_exon_peak*smooth_fraction)) or (f < cutoff_strict):
                    f_smooth = 0
                    cutoff = np.max([f, cutoff_, cutoff_strict])
                    last_pos = None
                else:
                    f_smooth = 1
                
                
        
        factor_smooth.append(f_smooth)
        
    last_pos = None
    for idx, f in reversed(list(enumerate(factor))):
        
        if (idx == (len(factor) - 1)):
            if (factor_smooth[idx]==1):
                last_pos = len(factor) - 1
            else:
                continue
        
        else:
            
            previous_f = factor_smooth[idx+1]
            current_f = factor_smooth[idx]
            
            if (previous_f == 0) and (current_f == 1):
                last_pos = idx
                continue
            
            elif (previous_f == 1):
                last_pos = np.min([last_pos, (idx+100)])
                
                if (current_f == 0) and (last_pos - idx) < 10:
                    # Remove small peaks that can be from mapping errors
                    for j in range(idx, last_pos+1):
                        factor_smooth[j] = 0
                    continue
                
                f_exon_peak = np.quantile(factor[(idx+1):(last_pos+1)], exon_quant)
                
                if (f > (f_exon_peak*smooth_fraction)) and (f > cutoff_strict):
                    factor_smooth[idx] = 1
                else:
                    continue

    factor_smooth = np.array(factor_smooth)
    
    return factor_smooth


def get_factors_stm(gene, k):
    
    EF = read_EF(gene, k)
    EF_binary = binarize_EF(EF)
    EF_shift = get_factor_shift(EF_binary)
    
    corrected_EF = pd.DataFrame()

    
    
    for i in range(1, k+1):
        factor = f'factor{str(i)}'
        factor_junctions_bed = make_factor_junctions_bed(EF_shift, factor=factor)
        junctions_match_list = get_factor_junctions(factor_junctions_bed)
        bed_factors_df, _ = get_factor_bed(EF_shift, factor)
        corrected_factor, idx_covered = correct_factor(
            EF, EF_binary, EF_shift, factor, junctions_match_list, bed_factors_df)
        
        corrected_EF[factor] = corrected_factor
        
    corrected_EF.index = EF.index
    
#     factors_bed = pd.DataFrame()
    
#     new_shift = get_factor_shift(corrected_EF)
#     for i in range(1, k+1):
        
    
    return EF, corrected_EF
        
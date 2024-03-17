import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import os
import subprocess as sp
from tqdm import tqdm 
import gtfparse

import seaborn as sns
import sys

sys.path.append('/project2/mstephens/cfbuenabadn/gtex-stm/code/scripts')
from plot_gtf import *

def plot_factors(EF, K, gene, gtf_df, plot_type='standard'):
    
    gtf_exons, cds = get_gene_gtf(gtf_df, gene)
    gene_name = gtf_df.loc[gtf_df.gene_id==gene].gene_name.iloc[0]
    
    if plot_type == 'standard':
        plot_gene(EF, gtf_df, gene, K, f'{gene_name}, fasttopics ini')
    
    elif plot_type == 'difference_track':
        plot_gene_difference(EF, gtf_df, gene, K, f'{gene_name}, fasttopics ini', difference=difference_track)
        
    elif plot_type == 'onerow':
        plot_gene_onerow(EF, gtf, gene, K, title= f'{gene_name}, fasttopics ini')
    
    
def plot_gene_onerow(EF, gtf, gene, K, title=None):
    gtf_exons, cds = get_gene_gtf(gtf_df, gene)
    colores = sns.color_palette("tab10")
    
    fig, ax = plt.subplots(nrows = 2, figsize=(10, 2*0.75), gridspec_kw={'height_ratios': [2, 3], 
                                                                           'wspace': 0.3, 'hspace': 0.3})

    start = int(EF.index[0].split(':')[1])
    end = int(EF.index[-1].split(':')[1])
    length = EF.shape[0]
    
    coords = np.linspace(start, end, num=length)
    
    for i in range(K):
        factor = f'factor{i+1}'
        scaled_y = EF[factor]/np.max(EF[factor])
        ax[0].plot(coords, scaled_y, c=colores[i], alpha=0.6)
        ax[0].set_xticks([])
        ax[0].spines['bottom'].set_visible(False)
        ax[0].spines['top'].set_visible(False)
        ax[0].spines['right'].set_visible(False)
        
        #ax[i].text(start+(length*0.9), 0.9, factor, size=12)
        
    if title:
        fig.suptitle(title, fontsize=12, x=0.5, y=1.1)

    PlotGene(gtf_exons, plot_cds=False, collapse_transcripts=False, plot_nmd=True, ax=ax[1])
    ax[1].set_xlim(ax[0].get_xlim())
    

def plot_gene_difference(EF, gtf, gene, K, title=None, difference = False):
    gtf_exons, cds = get_gene_gtf(gtf_df, gene)
    colores = sns.color_palette("tab10")
    
    if difference:
        nrows = K+2
    else:
        nrows = K+1
    
#     if K == 2:
#         S = 4*0.75
    if K <= 3:
        S = K
    else:
        S = K*0.8
    
    fig, ax = plt.subplots(nrows = nrows, figsize=(10, S), gridspec_kw={'height_ratios': [2]*(nrows-1) + [3], 
                                                                           'wspace': 0.3, 'hspace': 0.3})

    start = int(EF.index[0].split(':')[1])
    end = int(EF.index[-1].split(':')[1])
    length = EF.shape[0]
    
    coords = np.linspace(start, end, num=length)
    
    for i in range(K):
        factor = f'factor{i+1}'
        scaled_y = EF[factor]/np.max(EF[factor])
        ax[i].fill_between(coords,[0]*len(scaled_y), scaled_y, color=colores[i], alpha=0.9)
        ax[i].set_xticks([])
        ax[i].spines['bottom'].set_visible(False)
        ax[i].spines['top'].set_visible(False)
        ax[i].spines['right'].set_visible(False)
        
        ax[i].text(start+((end-start)*0.9), 0.9, factor, size=12)
        
    i = 2
    if difference:
        factor = 'difference'
        scaled_y = np.abs(np.array(EF.factor1) - np.array(EF.factor2))
        ax[i].fill_between(coords, [0]*len(scaled_y), scaled_y, color=colores[i], alpha=0.9)
        ax[i].set_xticks([])
        ax[i].spines['bottom'].set_visible(False)
        ax[i].spines['top'].set_visible(False)
        ax[i].spines['right'].set_visible(False)
        
        ax[i].text(start+((end-start)*0.9), 0.9, factor, size=12)
        
    if title:
        fig.suptitle(title, fontsize=12, x=0.5, y=1.1)

    PlotGene(gtf_exons, plot_cds=False, collapse_transcripts=False, plot_nmd=True, ax=ax[nrows-1])
    ax[nrows-1].set_xlim(ax[0].get_xlim())
    

def plot_coverage(gene, gtf, counts, tissue_list, title=None, log=True):
    gtf_exons, cds = get_gene_gtf(gtf_df, gene)
    colores = sns.color_palette("tab10")
        
    fig, ax = plt.subplots(nrows = (len(tissue_list)+1), figsize=(10, len(tissue_list)), 
                           gridspec_kw={'height_ratios': [2]*len(tissue_list) + [3], 
                                                                           'wspace': 0.3, 'hspace': 0.3})
    if log:
        counts = np.log1p(counts)

    start = int(counts.columns[0].split(':')[1])
    end = int(counts.columns[-1].split(':')[1])
    length = counts.shape[1]

    coords = np.linspace(start, end, num=length)
    
    for i, tissue in enumerate(tissue_list):
        
        tissue_idx = samples.loc[samples.tissue_id==tissue].index.intersection(counts.index)
    
        for idx, row in counts.loc[tissue_idx].iterrows():
            y = np.array(row)
            
            scaled_y = y/np.max(y)
            ax[i].plot(coords, scaled_y, c=colores[i], alpha=0.2, rasterized=True)
        
        ax[i].set_xticks([])
        ax[i].spines['bottom'].set_visible(False)
        ax[i].spines['top'].set_visible(False)
        ax[i].spines['right'].set_visible(False)
        
        if log:
            ax[i].set_ylabel('log1p\ncounts', size=9)
        else:
            ax[i].set_ylabel('counts', size=9)
        
        ax[i].text(start+((end-start)*0.65), 0.9, tissue, size=12)
        
    if title:
        fig.suptitle(title, fontsize=12)
    

    PlotGene(gtf_exons, plot_cds=False, collapse_transcripts=False, plot_nmd=True, ax=ax[len(tissue_list)])
    ax[len(tissue_list)].set_xlim(ax[0].get_xlim())
    
    
def plot_gene(EF, gtf_df, gene, K, title=None):
    gtf_exons, cds = get_gene_gtf(gtf_df, gene)
    colores = sns.color_palette("pastel")
    
    colores_diff = sns.color_palette("tab10")
    
    nrows = K+1
    
    if K <= 3:
        S = K
    else:
        S = K*0.8
    
    fig, ax = plt.subplots(nrows = nrows, figsize=(10, S), gridspec_kw={'height_ratios': [2]*(nrows-1) + [3], 
                                                                           'wspace': 0.3, 'hspace': 0.3})

    start = int(EF.index[0].split(':')[1])
    end = int(EF.index[-1].split(':')[1])
    length = EF.shape[0]
    
    coords = np.linspace(start, end, num=length)
    
    for i in range(K):
        factor = f'factor{i+1}'
        scaled_y = EF[factor]/np.max(EF[factor])
        ax[i].fill_between(coords,[0]*len(scaled_y), scaled_y, color=colores[i], alpha=0.9)
        
        list_other_factors = [f'factor{j+1}' for j in range(K) if (j != i)]
        
        sum_other_factors = np.zeros(len(scaled_y))
        
        for f in list_other_factors:
            sum_other_factors += np.array(EF[f])
            
        difference_y = scaled_y - sum_other_factors
        
        difference_y = np.maximum(difference_y, 0)
        
        common_y = scaled_y - difference_y
        
        ax[i].fill_between(coords,[0]*len(common_y), common_y, color=colores_diff[i], alpha=0.9)
        
        
        ax[i].set_xticks([])
        ax[i].spines['bottom'].set_visible(False)
        ax[i].spines['top'].set_visible(False)
        ax[i].spines['right'].set_visible(False)
        
        ax[i].text(start+((end-start)*0.9), 0.9, factor, size=12)
        
    
    if title:
        fig.suptitle(title, fontsize=12, x=0.5, y=1.1)

    PlotGene(gtf_exons, plot_cds=False, collapse_transcripts=False, plot_nmd=True, ax=ax[nrows-1])
    ax[nrows-1].set_xlim(ax[0].get_xlim())
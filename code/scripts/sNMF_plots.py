#!/usr/bin/env python

import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd

import sys
sys.path.append('/project2/mstephens/cfbuenabadn/gtex-stm/code')
sys.path.append('/project2/mstephens/cfbuenabadn/gtex-stm/code/scripts')
from get_unique_region_counts import *
from get_isoforms import *
from get_unique_region_counts import *
from collect_snmf_isoforms import *

'''
This script plots gtf genes and/or transcripts with matplotlib.
Unlike pyGenomeTracks, this script can be imported as a module 
to your python session to facilitate dynamic plotting.

There is still room for improvement, and I hope to keep adding to it
over time.

author: Carlos F Buen Abad Najar
email: cnajar@uchicago.edu
date: 04/11/2023
'''

def SmoothRectangle(x, y, width, height, x_smooth=0.1, smoothness_height=0.1, color='navy', ax=None):
    """
    A function to draw a smooth rectangle using Matplotlib.

    Args:
        x (float): The x-coordinate of the bottom left corner of the rectangle.
        y (float): The y-coordinate of the bottom left corner of the rectangle.
        width (float): The width of the rectangle.
        height (float): The height of the rectangle.
        x_smooth (float, optional): The amount of smoothing to apply to the x-axis. Defaults to 0.1.
        smoothness_height (float, optional): The amount of smoothing to apply to the y-axis. Defaults to 0.1.
        color (str, optional): The color of the rectangle. Defaults to 'navy'.
        ax (matplotlib.axes.Axes, optional): The axes on which to draw the rectangle. If not provided, a new figure and axes will be created. Defaults to None.

    Returns:
        matplotlib.patches.PathPatch: A PathPatch object representing the smooth rectangle.
    """
    if width*0.45 < x_smooth:
        x_smooth = 0

    y_smooth = height * smoothness_height
    
    verts = [
        (x + x_smooth, y),
        (x + width - x_smooth, y),
        (x + width, y),
        (x + width, y + y_smooth),
        (x + width, y + height - y_smooth),
        (x + width, y + height),
        (x + width - x_smooth, y + height),
        (x + x_smooth, y + height),
        (x, y + height),
        (x, y + height - y_smooth),
        (x, y + y_smooth),
        (x, y),
        (x + x_smooth, y)
    ]

    # Create a Path object from the vertex coordinates and codes
    Path = mpath.Path
    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.CURVE3,
        Path.CURVE3,
        Path.LINETO,
        Path.CURVE3,
        Path.CURVE3,
        Path.LINETO,
        Path.CURVE3,
        Path.CURVE3,
        Path.LINETO,
        Path.CURVE3,
        Path.CURVE3,
    ]
    path = Path(verts, codes)

    # If no axes are provided, create a new figure and axes
    if not ax:
        fig, ax = plt.subplots()
        
    # Create a PathPatch object from the path and add it to the axes
    patch = mpatches.PathPatch(path, facecolor=color, lw=2, edgecolor='none')
    ax.add_patch(patch)
    
    return patch


def PlotBedCoordinates(coords, y_position=0, height=1, smoothness=0.1, color='navy', ax=None):
    x, y = coords
    width = y - x
    if not ax:
        fig, ax = plt.subplots()
    SmoothRectangle(x, y_position, width, height, x_smooth=smoothness, color=color, ax=ax)


def PlotBed(coord_list, y_position=0, height=1, color='navy', ax=None):
    if not ax:
        fig, ax = plt.subplots()

    min_x = np.min([x for x, y in coord_list])
    max_x = np.max([y for x, y in coord_list])

    total_length = max_x - min_x

    alpha = total_length * 0.05

    smoothness = total_length * 0.005

    for coords in coord_list:
        fraction_plot = (coords[1] - coords[0]) / total_length
        PlotBedCoordinates(coords, y_position=y_position, height=height, smoothness=smoothness, color=color, ax=ax)


def PlotTranscriptCoords(transcript_gtf, position=0, cds_coords=None, color='navy', ax=None):
    # Extract coordinates from transcript GTF
    coord_list = [x for x in zip(transcript_gtf.start, transcript_gtf.end)]
    
    # Calculate the minimum and maximum x-coordinates
    min_x = np.min([x for x, y in coord_list])
    max_x = np.max([y for x, y in coord_list])
    
    # Calculate the total length of the transcript
    total_length = max_x - min_x
    
    # Set the alpha value
    alpha = total_length * 0.05
    
    # Create a new figure and axis if none are provided
    if not ax:
        fig, ax = plt.subplots(figsize=(15, 1))
        # Set the y-limits of the axis
        ax.set_ylim(-2.1, 3.1)
        # Set the x-limits of the axis
        ax.set_xlim(min_x - alpha, max_x + alpha)

    #############
    line_x = [np.min([int(x[1]) for x in coord_list])-1, np.max([int(x[0]) for x in coord_list])+1]
    ax.plot(line_x, [position + 0.5, position + 0.5], linewidth=1, color='black', zorder=0)
    ############################################
    # for arrows, uncomment this
    ############################################
    # Create blocks to represent the transcript

    ##############################
    ##############################
    
    # Plot the transcript bed
    PlotBed(coord_list, y_position=position, color=color, ax=ax)
    
    # Create a bed for the CDS coordinates if provided
    cds_bed = []
    
    if cds_coords:
        cds_x, cds_y = cds_coords
        cds_blocks = []
        
        # Iterate through the transcript coordinates and identify the CDS coordinates
        for x, y in coord_list:
            if y < cds_x:
                continue
            elif x > cds_y:
                break
            else:
                if cds_x > x:
                    x_coord = cds_x
                else:
                    x_coord = x
                if cds_y < y:
                    y_coord = cds_y
                else:
                    y_coord = y
                    
                coord_to_add = (x_coord, y_coord)
                cds_bed.append(coord_to_add)
                
        # Plot the CDS bed
        PlotBed(cds_bed, y_position=position - 0.5, height=2, color=color, ax=ax)


        
def PlotTranscript(gtf_transcript, position=0, cds=False, color='navy', ax=None):
    """
    Plots a transcript represented by a GTF file.

    Args:
        gtf_transcript (DataFrame): A DataFrame containing the GTF data for a transcript.
        position (int): The y-axis position of the transcript.
        cds (bool): Whether to highlight the coding sequence (CDS) of the transcript.
        color (str): The color to use for plotting the transcript.
        ax (matplotlib AxesSubplot): The AxesSubplot to use for plotting. If None, a new one is created.
    """
    # Create a new subplot if none is provided
    if not ax:
        fig, ax = plt.subplots(figsize=(15, 1))
        ax.set_ylim(-3.1, 3.1)
        
    # Determine the CDS coordinates if necessary
    if cds:
        try:
            try:
                start_codon = int(gtf_transcript.loc[gtf_transcript.feature == 'start_codon'].start)
            except:
                start_codon = int(gtf_transcript.loc[gtf_transcript.feature == 'start_codon'].start.min())

            try:
                stop_codon = int(gtf_transcript.loc[gtf_transcript.feature == 'stop_codon'].end)
            except:
                stop_codon = int(gtf_transcript.loc[gtf_transcript.feature == 'stop_codon'].end.max())

            cds = (start_codon, stop_codon)
        except:
            cds=None

    # Select the exons from the GTF file
    gtf_exons = gtf_transcript.loc[gtf_transcript.feature == 'exon']
    
    # Call the PlotTranscriptCoords function to plot the exons and CDS (if applicable)
    PlotTranscriptCoords(gtf_exons[['start', 'end', 'strand']], position=position, cds_coords=cds, color=color, ax=ax)


def PlotGTF_ax(gene_gtf, ax, plot_nmd=False, plot_cds=False, collapse_transcripts=True):
    # Filter out non-protein coding transcripts
    pc_gtf = gene_gtf.loc[gene_gtf.transcript_type == 'protein_coding']
    
    # Separate out NMD and retained intron transcripts
    nmd_gtf = gene_gtf.loc[gene_gtf.transcript_type.isin(['nonsense_mediated_decay', 'retained_intron'])]
    
    # Get unique protein coding transcripts
    pc_transcripts = pc_gtf.transcript_name.unique()
    
    # Initialize position
    position = 0
    
    # Loop over protein coding transcripts
    for transcript, transcript_gtf in pc_gtf.groupby('transcript_name'):
        # Plot transcript
        PlotTranscript(transcript_gtf, position = position, ax=ax, cds=plot_cds)
        
        # Adjust position for next transcript
        if not collapse_transcripts:
            if plot_cds:
                position -= 2.5
            else:
                position -= 1.5
                
    # If plot_nmd is True, plot NMD transcripts
    if plot_nmd:
        # Set color for NMD transcripts depending on whether or not transcripts are being collapsed
        if collapse_transcripts:
            color = 'navy'
        else:
            color = 'darkorange'
            
        # Loop over NMD transcripts
        for transcript, transcript_gtf in nmd_gtf.groupby('transcript_name'):
            # Plot transcript
            PlotTranscript(transcript_gtf, position = position, ax=ax, cds=False, color=color)

            # Adjust position for next transcript
            if not collapse_transcripts:
                if plot_cds:
                    position -= 2.5
                else:
                    position -= 1.5
                
    # Set y-axis limits
    ax.set_ylim(position, 2)

    # Remove spines and ticks
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])


def PlotGene(gene_gtf, plot_nmd=False, plot_cds=False, collapse_transcripts=True, ax=None):
    pc_gtf = gene_gtf.loc[gene_gtf.transcript_type == 'protein_coding']
    nmd_gtf = gene_gtf.loc[gene_gtf.transcript_type.isin(['nonsense_mediated_decay', 'retained_intron'])]
    
    pc_transcripts = pc_gtf.transcript_name.unique()
    
    if not ax:
        fig, ax = plt.subplots(figsize=(15, 1))
        ax.set_ylim(-6.5, 2.1)
    
    position = 0
    
    for transcript, transcript_gtf in pc_gtf.groupby('transcript_name'):
        PlotTranscript(transcript_gtf, position = position, ax=ax, cds=plot_cds)
        
        if not collapse_transcripts:
            if plot_cds:
                position -= 2.5
            else:
                position -= 1.5
                
    if plot_nmd:
        if collapse_transcripts:
            color = 'navy'
        else:
            color = 'navy'
            
        for transcript, transcript_gtf in nmd_gtf.groupby('transcript_name'):
            PlotTranscript(transcript_gtf, position = position, ax=ax, cds=False, color=color)

            if not collapse_transcripts:
                if plot_cds:
                    position -= 2.5
                else:
                    position -= 1.5
                
    

    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_linewidth(1)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    ticks_loc = ax.get_xticks().tolist()
    label_format = '{:,}'

    ax.set_xticks(ticks_loc)
    ax.set_xticklabels([label_format.format(int(x)) for x in ticks_loc])
    
    chrom = list(gene_gtf.seqname)[0]
    
    ax.set_xlabel(chrom + ' coordinates')
    ax.set_yticks([])
        
def get_gene_gtf(gtf_df, gene):
    try:
        gtf_gene = gtf_df.loc[gtf_df.gene_id==gene]
    except:
        gtf_gene = gtf_df.loc[gtf_df.gene_name==gene]
    gtf_exons = gtf_gene.loc[gtf_gene.feature=='exon']
    strand = gtf_gene.strand.iloc[0]
    
    if strand == '+':
        start = gtf_gene.loc[gtf_gene.feature=='start_codon'].start.min()
        stop = gtf_gene.loc[gtf_gene.feature=='stop_codon'].end.max()
    else:
        start = gtf_gene.loc[gtf_gene.feature=='stop_codon'].start.min()
        stop = gtf_gene.loc[gtf_gene.feature=='start_codon'].end.max()
    
    cds = (start, stop)
    
    return gtf_exons, cds

def plot_factor_tracks(EF, gene, K, title=None, fill=True, smooth=False, colores = None, q=0.99, figsize=None):
    # gtf_exons, cds = get_gene_gtf(gtf_df, gene)

    EF = np.minimum(EF/EF.quantile(q, axis=0), 1)

    if colores is None:
        if K <= 5:
            colores = ['tab:blue', 'tab:orange', 'tab:green', 'goldenrod', 'tab:red'] 
        else:
            colores = sns.color_palette("tab10")
    
    if K <= 3:
        S = K
    else:
        S = K*0.8

    if figsize is None:
        figsize=(15, S)
    
    fig, ax = plt.subplots(nrows = K, figsize=figsize, gridspec_kw={'height_ratios': [2]*K, 
                                                                           'wspace': 0.3, 'hspace': 0.3})

    start = int(EF.index[0].split(':')[1])
    end = int(EF.index[-1].split(':')[1])
    length = EF.shape[0]
    
    coords = np.linspace(start, end, num=length)
    
    for i in range(K):
        factor = f'factor{i+1}'
        scaled_y = EF[factor]#/np.max(EF[factor])
        if fill:
            ax[i].fill_between(coords, np.zeros(len(coords)), scaled_y, color=colores[i], alpha=0.3)
#         else:
        ax[i].plot(coords, scaled_y, c=colores[i], alpha=0.9, linewidth=3)
        ax[i].set_xticks([])
        ax[i].spines['bottom'].set_visible(False)
        ax[i].spines['top'].set_visible(False)
        ax[i].spines['right'].set_visible(False)
        
        ax[i].text(start+((end-start)*1.1), 0.9, factor, size=12)
        
    if title:
        fig.suptitle(title, fontsize=12, x=0.5, y=1.1)

    # PlotGene(gtf_exons, plot_cds=False, collapse_transcripts=False, plot_nmd=True, ax=ax[K])
    # ax[K].set_xlim(ax[0].get_xlim())
    


def plot_gene_(EF, gtf_df, gene, K, title=None, fill=True, smooth=False, colores = None):
    gtf_exons, cds = get_gene_gtf(gtf_df, gene)

    if colores is None:
        if K <= 5:
            colores = ['tab:blue', 'tab:orange', 'tab:green', 'goldenrod', 'tab:red'] 
        else:
            colores = sns.color_palette("tab10")
    
    if K <= 3:
        S = K
    else:
        S = K*0.8
    
    fig, ax = plt.subplots(nrows = K+1, figsize=(10, S), gridspec_kw={'height_ratios': [2]*K + [3], 
                                                                           'wspace': 0.3, 'hspace': 0.3})

    start = int(EF.index[0].split(':')[1])
    end = int(EF.index[-1].split(':')[1])
    length = EF.shape[0]
    
    coords = np.linspace(start, end, num=length)
    
    for i in range(K):
        factor = f'factor{i+1}'
        scaled_y = EF[factor]/np.max(EF[factor])
        if fill:
            ax[i].fill_between(coords, np.zeros(len(coords)), scaled_y, color=colores[i], alpha=0.3)
#         else:
        ax[i].plot(coords, scaled_y, c=colores[i], alpha=0.9, linewidth=2)
        ax[i].set_xticks([])
        ax[i].spines['bottom'].set_visible(False)
        ax[i].spines['top'].set_visible(False)
        ax[i].spines['right'].set_visible(False)
        
        ax[i].text(start+((end-start)*1.1), 0.9, factor, size=12)
        
    if title:
        fig.suptitle(title, fontsize=12, x=0.5, y=1.1)

    PlotGene(gtf_exons, plot_cds=False, collapse_transcripts=False, plot_nmd=True, ax=ax[K])
    ax[K].set_xlim(ax[0].get_xlim())
    

def plot_coverage(gene, gtf_df, counts, tissue_list, title=None, log=True):
    gtf_exons, cds = get_gene_gtf(gtf_df, gene)
    colores = sns.color_palette("tab10")
    
#     tissues = sorted(os.listdir('coverage/counts_tables/'))
    
    fig, ax = plt.subplots(nrows = (len(tissue_list)+1), figsize=(20, 2*len(tissue_list)), 
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
        
        ax[i].text(start+((end-start)*0.95), 0.9, tissue, size=12)
        
    if title:
        fig.suptitle(title, fontsize=12)
    

    PlotGene(gtf_exons, plot_cds=False, collapse_transcripts=False, plot_nmd=True, ax=ax[len(tissue_list)])
    ax[len(tissue_list)].set_xlim(ax[0].get_xlim())
    
#     plt.show()

def plot_coverage(gene, gtf_df, counts, tissue, color, title=None, log=True, height=2.5):
    gtf_exons, cds = get_gene_gtf(gtf_df, gene)
    
    
#     tissues = sorted(os.listdir('coverage/counts_tables/'))
    
    fig, ax = plt.subplots(nrows = 2, figsize=(10, height), 
                           gridspec_kw={'height_ratios': [2, 3], 'wspace': 0.3, 'hspace': 0.3})
    if log:
        counts = np.log1p(counts)

    start = int(counts.columns[0].split(':')[1])
    end = int(counts.columns[-1].split(':')[1])
    length = counts.shape[1]

    coords = np.linspace(start, end, num=length)
    
#     for i, tissue in enumerate(tissue_list):
        
    tissue_idx = samples.loc[samples.tissue_id==tissue].index.intersection(counts.index)

#     tissue_quant = counts.loc[tissue_idx].max(axis=1).quantile(0.99)
    
    scaled_y = np.array(counts.loc[tissue_idx].max(axis=0))
    
    tissue_quant = np.max(scaled_y)* 0.95

    for idx, row in counts.loc[tissue_idx].iterrows():
        y = np.array(row)

        scaled_y = y#/np.max(y)
        ax[0].plot(coords, scaled_y, c=color, alpha=0.5, linewidth=0.5, rasterized=True)

    ax[0].set_xticks([])
    ax[0].spines['bottom'].set_visible(False)
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['right'].set_visible(False)

    if log:
        ax[0].set_ylabel('log1p\ncounts', size=9)
    else:
        ax[0].set_ylabel('counts', size=9)

    ax[0].text(start+((end-start))*1.1, 0.7*(np.max(scaled_y)), tissue, size=12, c = color)

    ax[0].set_ylim([-0.1, tissue_quant])
        
    if title:
        fig.suptitle(title, fontsize=12)
    

    PlotGene(gtf_exons, plot_cds=False, collapse_transcripts=False, plot_nmd=True, ax=ax[1])
    ax[1].set_xlim(ax[0].get_xlim())

def plot_gene(EF, gtf_df, gene, K, title=None):
    gtf_exons, cds = get_gene_gtf(gtf_df, gene)
    colores = ['tab:blue',  'tab:orange', 'tab:green', 'goldenrod', 'tab:red']#sns.color_palette("pastel")
    
    colores_diff = colores#['tab:blue', 'tab:orange',  'tab:green', 'goldenrod', 'tab:red']#sns.color_palette("tab10")
    
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
        ax[i].fill_between(coords,[0]*len(scaled_y), scaled_y, color=colores[i], alpha=0.5)
        
        list_other_factors = [f'factor{j+1}' for j in range(K) if ((j != i) and (j!=1))]
        
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
        
        ax[i].text(start+((end-start)*1.1), 0.9, factor, size=12)
        
    
    if title:
        fig.suptitle(title, fontsize=12, x=0.5, y=1.1)

    PlotGene(gtf_exons, plot_cds=False, collapse_transcripts=False, plot_nmd=True, ax=ax[nrows-1])
    ax[nrows-1].set_xlim(ax[0].get_xlim())

def plot_factors_wrapper(rds, gene, ebpmf='ebpmf_10', strand='plus', colores=None, figsize=None, lm=True, factor_list = None):
    EF = pd.DataFrame(rds[ebpmf]['train_fit']['EF_smooth'][1:-1,:])
    EF.index = rds[ebpmf]['coords']
    EF.columns = [f'factor{str(i+1)}' for i in range(EF.shape[1])]

    if factor_list is not None:
        EF = EF[factor_list]
        EF = pd.DataFrame(EF)
        K = EF.shape[1]
        EF.columns = [f'factor{str(i)}' for i in range(1, K+1)]
    
    for factor in EF.columns:
        
        if lm:
            y = factor_lm(EF[factor], strand)
        else:
            y = np.array(EF[factor])
        EF[factor] = y
    plot_factor_tracks(EF/EF.quantile(0.99, axis=0), gene, EF.shape[1], colores = colores, figsize=figsize)

def plot_isoform_annotations(annotation_exons, gene, colores=None, start=None, end=None, figsize=None, lwidth=5, iso_order=None):

    gene_exons = annotation_exons.loc[annotation_exons.gene_id == gene]

    if iso_order is None:
        isoforms = sorted(gene_exons.transcript_id.unique())
    else:
        isoforms = [gene + '.' + x for x in iso_order]

    isoform_dict = {}
    for i, iso in enumerate(isoforms):
        isoform_name = f'isoform_{str(i+1)}'
        df = gene_exons.loc[gene_exons.transcript_id == iso]
        df['transcript_id'] =  f'{gene}.{isoform_name}'
        isoform_dict.update({isoform_name:{'df':df}})

    chrom = list(annotation_exons.chrom)[0]
    if start is None:
        start = str(np.min([int(list(gene_exons.start)[0]), int(list(gene_exons.start)[0])]) - 1000)
    if end is None:
        end = str(np.max([int(list(gene_exons.end)[-1]), int(list(gene_exons.end)[-1])]) + 1000)

    coords = [f'{chrom}:{start}', f'{chrom}:{end}']

    # print(isoform_dict)

    plot_gene_isoforms(isoform_dict, coords, color_list = colores, figsize=figsize, lwidth=lwidth)


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


def plot_gene_isoforms(isoforms_dict, coordinates, color_list = None, axes=None, figsize=None, lwidth=5):
    xlim1 = int(coordinates[0].split(':')[1])
    xlim2 = int(coordinates[-1].split(':')[1])

    if color_list is None:
        color_list = sns.color_palette("tab10")

    K = len(isoforms_dict)

    if figsize is None:
        figsize=(20, 3)

    if axes is None:
        fig, axes = plt.subplots(K, 1, figsize=figsize)

    for i in range(K):
        isoform_df = isoforms_dict[f'isoform_{str(i+1)}']['df']
        # print(isoform_df)
        ax = axes[i]
        color = color_list[i]
        plot_isoform(isoform_df, ax, color, xlim1, xlim2, lwidth=lwidth)
        

def plot_isoform(isoform_df, ax, color, xlim1, xlim2, lwidth):
    is_first = True
    for idx, row in isoform_df.iterrows():
        start = int(row.start)
        end = int(row.end)
        if is_first:
            first = end
            is_first = False
    
        ax.fill_between([start, end], [0, 0], [1, 1], color = color, zorder=2)
    
    ax.plot([first, start], [0.5, 0.5], c=color, linewidth=lwidth)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines[['bottom', 'top', 'right', 'left']].set_visible(False)
    ax.set_xlim([xlim1, xlim2])



def plot_bar_ax(loadings, color, ax, label='label', sort_factor = False):
    x = range(len(loadings))
    if sort_factor:
        loadings = sorted(loadings) 
    ax.bar(x, loadings, width=4, color=color, alpha=0.7)
    if sort_factor:
        median = np.median(loadings)
        middle_point = len(loadings)/2
        ax.plot([0, len(loadings)], [median, median], linestyle = '--', c='black', linewidth=2)
        ax.scatter([middle_point], [median], marker="D", c = 'black', s=10)
    
    
    
def factor_barplot(EL, samples, color_list, label_list, sort_factor = False, figsize=None, ylim = None):
    samples = samples.loc[EL.index]
    tissues = samples.tissue_id.unique()
    factors = EL.columns
    
    K = len(factors)
    N = len(tissues)

    if figsize is None:
        figsize = (N/3,K/1.2)
    
    
    

    if K == 2:
        fig, axes = plt.subplots(1,N+1, figsize=figsize, gridspec_kw={'wspace': 0.1, 'hspace': 0.1}, 
                             width_ratios = ([1]+([4]*N)))

        for j, tissue in enumerate(tissues):
            label = label_list[j]
            tissue_samples = samples.loc[samples.tissue_id == tissue].index
            loadings = list(EL.loc[tissue_samples, factors[0]])

            ax = axes[j+1]

            loadings = sorted(loadings)

            median = np.median(loadings)
            middle_point = len(loadings)/2
            for i, loading in enumerate(loadings):
                ax.plot([i, i], [0, loading], c=color_list[0])
                ax.plot([i, i], [loading, 1], c=color_list[1])
                
            ax.plot([0, len(loadings)], [median, median], linestyle = '--', c='black', linewidth=2)
            ax.scatter([middle_point], [median], marker="D", c = 'black', s=10, zorder=len(loadings)+2)

            axes[0].spines['top'].set_visible(False)
            axes[0].spines['right'].set_visible(False)
            axes[0].spines['bottom'].set_visible(False)
            axes[0].set_xticks([])
            ax.set_xlabel(label, rotation=45)

            if ylim is None:
                axes[0].set_ylim([0,1])

                ax.set_ylim([0,1])
            else:
                axes[0].set_ylim(ylim)
                ax.set_ylim(ylim)
            ax.spines['top'].set_visible(False)
            
            ax.spines['left'].set_visible(False)
            ax.set_yticks([])
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.set_xticks([])
            # ax[0].set_ylabel(factor)
            #plot_bar_ax(loadings, color, ax, label=label, sort_factor=True)

    else:
        fig, axes = plt.subplots(K,N+1, figsize=figsize, gridspec_kw={'wspace': 0.1, 'hspace': 0.1}, 
                             width_ratios = ([1]+([4]*N)))
    
        for i, factor in enumerate(factors):
            color = color_list[i]
            
            axes[i,0].spines['top'].set_visible(False)
            axes[i,0].spines['right'].set_visible(False)
            axes[i,0].spines['bottom'].set_visible(False)
            axes[i,0].set_xticks([])
            axes[i,0].set_ylabel(factor)
            
            factor_max = EL[factor].max()
            
            axes[i,0].set_ylim([0,factor_max])
            
            for j, tissue in enumerate(tissues):
                if i == (len(factors)-1):
                    label = label_list[j]
                else:
                    label = ''
                tissue_samples = samples.loc[samples.tissue_id == tissue].index
                loadings = list(EL.loc[tissue_samples, factor])
                
                    
                ax = axes[i,j+1]
                plot_bar_ax(loadings, color, ax, label=label, sort_factor=sort_factor)
                
                ax.set_ylim([0,factor_max])
                ax.spines['top'].set_visible(False)
                
                ax.spines['left'].set_visible(False)
                ax.set_yticks([])
                ax.spines['right'].set_visible(False)
                ax.spines['bottom'].set_visible(False)
                ax.set_xticks([])
                
    
            #     ax.axis('off')  # Hide
                ax.margins(0, 0)
                
    
                ax.set_xlabel(label, rotation=45)



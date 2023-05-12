import argparse

import pandas as pd
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import numpy as np
import gtfparse
import rpy2.robjects as robjects

from rpy2.robjects import pandas2ri
import seaborn as sns
import sys

sys.path.append('/project2/mstephens/cfbuenabadn/gtex-stm/code/scripts')
from plot_gtf import *

import rpy2.robjects.conversion

# Define a Python function to read RDS files
def read_rds(filename):
    with robjects.conversion.localconverter(robjects.default_converter + robjects.pandas2ri.converter):
        return robjects.r['readRDS'](filename)
    


def read_ebpmf_object(ebpmf_rds):
    ebpmf_object = read_rds(ebpmf_rds)

    EL = ebpmf_object['fit_ebpmf']['EL']
    K = EL.shape[1]

    samples = np.array(ebpmf_object['samples'])
    annotation = pd.DataFrame([x.split('.') for x in samples], columns=['tissue_id', 'sample_id']).set_index('sample_id')
    
    coords = np.array(ebpmf_object['coords'])
    counts = pd.DataFrame(ebpmf_object['geneCounts'], columns=coords, index=annotation.index).T

    # get the number of zeros in the counts dataframe and slice the coords array accordingly
    counts_zero = counts.sum(axis=1).values
    K_left = np.argmax(counts_zero > 0)
    K_right = np.argmax(counts_zero[::-1] > 0)
    
    coords_for_EF = coords[K_left:-K_right] if K_right > 0 else coords[K_left:]

    EL = pd.DataFrame(EL, index=annotation.index, columns=['factor'+str(i+1) for i in range(K)])
    EF = pd.DataFrame(ebpmf_object['fit_ebpmf']['EF'], index=coords_for_EF, columns=['factor'+str(i+1) for i in range(K)])
    LF = EF.dot(EL.T)

    counts['coords'] = [int(x.split('.')[1]) for x in counts.index]
    EF['coords'] = [int(x.split('.')[1]) for x in EF.index]
    LF['coords'] = [int(x.split('.')[1]) for x in LF.index]

    annotation = pd.DataFrame([x.split('.') for x in samples], columns=['tissue_id', 'sample_id']).set_index('sample_id')
    
    return counts, EF, EL, LF, annotation



def plot_line_tracks(coords, track, ax, **kwargs):
    
    coords = np.array(coords)
    track = np.array(track)
    
    if len(track.shape) == 1:
        ax.plot(coords, track, **kwargs)
    else:
        tracks_to_plot = track.shape[1]
        
        for i in range(tracks_to_plot):
            ax.plot(coords, track[:, i], alpha=0.05, linewidth=0.2, **kwargs)
            
        mean = np.mean(track, axis=1)
        ax.plot(coords, mean, linewidth=1, **kwargs)
        
    
def plot_factors(EF, gtf, gene_name='', rolling_mean = 43, savepdf = None):
    
    EF = EF.rolling(rolling_mean).mean()
    
    colores = sns.color_palette("Paired", 10)
    
    
    coords = np.array(EF['coords'])
    
    factors = [x for x in EF.columns if 'factor' in x]
    
    K = len(factors)
    
    transcript_total = gtf.shape[1]
    
    gtf_ax_ysize = np.max([int(transcript_total / 5), 1])
    
    height_ratios = ([1] * (K)) + [0.3, 0.1, gtf_ax_ysize]
    
    K_ax_ysize = K*(1.3)
    
    
    fig, ax = plt.subplots(nrows = (K + 3), figsize=(15, K_ax_ysize + gtf_ax_ysize), 
                           gridspec_kw={'height_ratios': height_ratios, 
                                                                       'wspace': 0.3, 'hspace': 0.6})

    PlotGTF_ax(gtf, plot_cds=False, plot_nmd=True, collapse_transcripts=False, ax=ax[-1])
    
    for i in range(K):
        factor_to_plot = np.array(EF[factors[i]])
        color = colores[i]
        plot_line_tracks(coords, factor_to_plot, ax=ax[i], c=colores[i], linewidth=2)
        
        
        ax[i].spines['right'].set_visible(False)
        ax[i].spines['bottom'].set_visible(False)
        ax[i].spines['top'].set_visible(False)
        ax[i].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax[i].text(0.01, 0.3, 'Factor ' + str(i+1), transform=ax[i].transAxes, size=14)
        if i < K-1:
            ax[i].set_xticks([])
            
        else:
            
            ticks_loc = ax[i].get_xticks().tolist()

            label_format = '{:,}'

            ax[i].set_xticks(ticks_loc)
            ax[i].set_xticklabels([label_format.format(int(x)) for x in ticks_loc])

            xlim = ax[i].get_xlim()
            
    ax[K-1].set_xticks([])
            
    ax[K].spines['left'].set_visible(False)
    ax[K].spines['bottom'].set_linewidth(2)
    ax[K].spines['top'].set_visible(False)
    ax[K].spines['right'].set_visible(False)
    ax[K].set_yticks([])
    ax[K].set_xticks(ticks_loc)
    ax[K].set_xticklabels([label_format.format(int(x)) for x in ticks_loc], size=14)
    chrom = list(gtf.seqname)[0]
    ax[K].set_xlabel(chrom + ' coordinates', size=14)
    ax[K].patch.set_visible(False)
    ax[K].xaxis.set_tick_params(width=2)
    
    ax[K+1].spines['left'].set_visible(False)
    ax[K+1].spines['bottom'].set_visible(False)
    ax[K+1].spines['top'].set_visible(False)
    ax[K+1].spines['right'].set_visible(False)
    ax[K+1].set_xticks([])
    ax[K+1].set_yticks([])

    ax[K+1].patch.set_visible(False)

    
    for ax_ in ax:
        ax_.set_xlim(xlim)
        
    ax[K+2].text(0.01, 0.9, gene_name, transform=ax[K+2].transAxes, size=14)
        
    if savepdf:
        plt.savefig(savepdf, 
                   bbox_inches='tight', 
                   transparent=True,
                   pad_inches=0)
        

def plot_samples_track(counts, annotation, gtf, gene_name = '', count_data=True, rolling_mean=43, savepdf = None):
    
    if count_data:
        coords = np.array(counts.coords)
        counts = np.log1p(counts)
        counts['coords'] = coords
    
    counts = counts.rolling(rolling_mean).mean()
    
    colores_tissue = plt.get_cmap('tab10').colors

    
    tissues = sorted(annotation.loc[counts.drop('coords', axis=1).columns].tissue_id.unique())
    
    coords = np.array(counts['coords'])
    
    K = len(tissues)
    
    transcript_total = gtf.shape[1]
    
    gtf_ax_ysize = np.max([int(transcript_total / 5), 1])
    
    height_ratios = ([1] * (K)) + [0.3, 0.3, gtf_ax_ysize]
    
    K_ax_ysize = K*(1.5)
    
    
    fig, ax = plt.subplots(nrows = (K + 3), figsize=(15, K_ax_ysize + gtf_ax_ysize), 
                           gridspec_kw={'height_ratios': height_ratios, 'wspace': 0.3, 'hspace': 0.6})
    
    PlotGTF_ax(gtf, plot_cds=False, plot_nmd=True, collapse_transcripts=False, ax=ax[-1])

    for i in range(K):
        tissue = tissues[i]
        
        tissue_counts = counts[annotation.loc[annotation.tissue_id == tissue].index]
        if count_data:
            tissue_counts = np.log1p(tissue_counts)
            
        plot_line_tracks(coords, tissue_counts, ax=ax[i], c=colores_tissue[i])
        
        ax[i].spines['right'].set_visible(False)
        ax[i].spines['bottom'].set_visible(False)
        ax[i].spines['top'].set_visible(False)
        ax[i].ticklabel_format(style='sci', axis='y', scilimits=(0,0))

        
        label = ' '.join(tissue.split('_'))
        if count_data:
            label += ' (counts)'
        else:
            label += r' ($LF$)'
        ax[i].text(0.1, 1.0, label, transform=ax[i].transAxes, size=14)
        if i < K-1:
            ax[i].set_xticks([])
            
        else:
            
            ticks_loc = ax[i].get_xticks().tolist()

            label_format = '{:,}'

            ax[i].set_xticks(ticks_loc)
            ax[i].set_xticklabels([label_format.format(int(x)) for x in ticks_loc])

            xlim = ax[i].get_xlim()
            
    ax[K-1].set_xticks([])
            
    ax[K].spines['left'].set_visible(False)
    ax[K].spines['bottom'].set_linewidth(2)
    ax[K].spines['top'].set_visible(False)
    ax[K].spines['right'].set_visible(False)
    ax[K].set_yticks([])
    ax[K].set_xticks(ticks_loc)
    ax[K].set_xticklabels([label_format.format(int(x)) for x in ticks_loc], size=14)
    chrom = list(gtf.seqname)[0]
    ax[K].set_xlabel(chrom + ' coordinates', size=14)
    ax[K].patch.set_visible(False)
    ax[K].xaxis.set_tick_params(width=2)
    
    ax[K+1].spines['left'].set_visible(False)
    ax[K+1].spines['bottom'].set_visible(False)
    ax[K+1].spines['top'].set_visible(False)
    ax[K+1].spines['right'].set_visible(False)
    ax[K+1].set_xticks([])
    ax[K+1].set_yticks([])

    ax[K+1].patch.set_visible(False)

    
    for ax_ in ax:
        ax_.set_xlim(xlim)
        
    ax[K+2].text(0.01, 0.9, gene_name, transform=ax[K+2].transAxes, size=14)
        
    if savepdf:

        plt.savefig(savepdf, 
                   bbox_inches='tight', 
                   transparent=True,
                   pad_inches=0)
    
parser = argparse.ArgumentParser()
parser.add_argument('--GTF', type=str, required=True)
parser.add_argument('--RDS', type=str, required=True)
parser.add_argument('--geneName', type=str, required=True)
parser.add_argument('--kfactors', type=str, required=True)
parser.add_argument('FILES', type=str, nargs="+")

if __name__ == '__main__':

    args = parser.parse_args()
    
    gtf_file = args.GTF
    rds_file = args.RDS
    geneName = args.geneName
    kfactors = args.kfactors
    output = args.FILES
    
    factors_plot, counts_plot, LF_plot = output
    
    # Parse GTF file to DataFrame
    gtf_df = gtfparse.read_gtf(gtf_file)

    gene_gtf = gtf_df.loc[(gtf_df.feature == 'exon') & (gtf_df.gene_name == geneName)]
    
    #rds_template = "/project2/mstephens/cfbuenabadn/gtex-stm/code/ebpmf_model/train_10tissues/models/{geneName}.K{kfactors}.ebpmf.rds"
    #rds_file = rds_template.format(geneName=geneName, kfactors=kfactors)
    
    counts, EF, EL, LF, annotation = read_ebpmf_object(rds_file)
    
    plot_factors(EF, gene_gtf, gene_name = geneName, savepdf = factors_plot)
    
    plot_samples_track(counts, annotation, gene_gtf, count_data=True, gene_name = geneName, savepdf = counts_plot)
    
    plot_samples_track(LF, annotation, gene_gtf, count_data=False, gene_name = geneName, savepdf = LF_plot)


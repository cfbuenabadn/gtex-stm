#!/usr/bin/env python

import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import numpy as np

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
    
    # Create blocks to represent the transcript
    blocks = np.linspace(np.min([int(x[1]) for x in coord_list])-1, np.max([int(x[0]) for x in coord_list])+1, num=101)

    # Calculate the step size
    step = blocks[1] - blocks[0]

    # Create arrows to represent the blocks
    for i in range(100):
        plt.arrow(blocks[i], position + 0.5, step, 0, head_width=0.25, head_length=step/3, fc='k', 
                  ec='k', linestyle='-', linewidth=1, color='grey',
                  alpha=0.3, zorder=0)

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
    PlotTranscriptCoords(gtf_exons[['start', 'end']], position=position, cds_coords=cds, color=color, ax=ax)


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

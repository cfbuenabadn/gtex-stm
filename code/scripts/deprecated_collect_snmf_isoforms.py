import numpy as np
import pandas as pd
import rpy2
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
import rpy2.robjects.packages as rpackages
ebpmf = rpackages.importr('ebpmf.alpha')
import os

from matplotlib import pyplot as plt
import seaborn as sns
import tabix
from pybedtools import BedTool
import sys

sys.path.append('/project2/mstephens/cfbuenabadn/gtex-stm/code/scripts')
from get_isoforms import *
from prepare_counts import *
from sNMF_plots import *
from tqdm import tqdm
import pingouin as pg

def load_ebpmf_gene(gene_id):
    readRDS = ro.r['readRDS']
    df = readRDS(f'ebpmf_models/filtered/RDS/{gene_id}.rds')
    with (ro.default_converter + pandas2ri.converter).context():
        pd_from_r_df = ro.conversion.get_conversion().rpy2py(df)

    junc_file = '/project2/mstephens/cfbuenabadn/gtex-stm/code/junctions.tab.gz'
    coords = pd_from_r_df['coords']
    junctions_bed = run_tabix_on_junc(junc_file, coords, gene_id)

    output = {'rds':pd_from_r_df, 'junctions_bed':junctions_bed}

    return output

genes = os.listdir('ebpmf_models/filtered/RDS/')
genes = [x.split('.')[0] for x in genes]
samples = pd.read_csv('config/samples.tsv', sep='\t', index_col=0)

def write_gtf(merged_isoforms, gene, strand):
    K = len(merged_isoforms)
    isoforms = [f'isoform_{str(i+1)}' for i in range(K)]
    gene_start = np.nan
    gene_end = np.nan
    chrom = merged_isoforms['isoform_1']['df'].chrom.iloc[0]

    gtf_lines = []

    gene_description = f'gene_id "{gene}"'
    
    for isoform in isoforms:
        isoform_name = gene + '.' + isoform
        df = merged_isoforms[isoform]['df']
        iso_start = df.start.astype(int).min()
        iso_end = df.end.astype(int).max()

        gene_start = np.nanmin([gene_start, iso_start])
        gene_end = np.nanmax([gene_end, iso_end])

        iso_start = str(iso_start)
        iso_end = str(iso_end)

        isoform_description = gene_description + f'; transcript_id "{isoform_name}"'


        isoform_row = [chrom, 's-NMF', 'transcript', iso_start, iso_end, '.', strand, '.', isoform_description]
        isoform_row = '\t'.join(isoform_row)

        gtf_lines.append(isoform_row)

        exon_count = 1
        for idx, row in df.iterrows():
            exon_start = str(row.start)
            exon_end = str(row.end)
            exon_name = f'exon_{str(exon_count)}'
            exon_description = isoform_description + f'; exon_id "{exon_name}"'
            exon_row = [chrom, 's-NMF', 'exon', exon_start, exon_end, '.', strand, '.', exon_description]
            exon_row = '\t'.join(exon_row)
            gtf_lines.append(exon_row)
            exon_count += 1

    gene_start = str(gene_start)
    gene_end = str(gene_end)

    gene_line = [chrom, 's-NMF', 'gene', gene_start, gene_end, '.', strand, '.', gene_description]
    gene_line = '\t'.join(gene_line)

    gtf_lines = [gene_line] + gtf_lines

    gene_gtf = '\n'.join(gtf_lines) + '\n'
    
    return gene_gtf 

def run_anova(EL_manova, isoform):

    res = pg.anova(dv=isoform, 
                 between='tissue_id', 
                 data=EL_manova,
                 detailed=True)

    SS_b = res.iloc[0, 1]
    SS_tot = res.iloc[1, 1]
    eta_squared = SS_b/SS_tot
    try:
        pvals = res.iloc[0, 5]
    except IndexError:
        return 'NA', 'NA', 'NA', 'NA'

    kw = pg.kruskal(dv=isoform, 
         between='tissue_id', 
         data=EL_manova,
         detailed=True)

    H = kw.iloc[0, 2]
    pvals_kw = kw.iloc[0, 3]

    return eta_squared, pvals, H, pvals_kw

def merge_EL(EL, merged_isoforms):
    K = EL.shape[1]
    EL = pd.DataFrame(EL)
    EL.columns = [f'factor_{str(i+1)}' for i in range(K)]
    EL_merged = pd.DataFrame()
    K_ = len(merged_isoforms)

    for i in range (K_):
        iso_name = f'isoform_{str(i+1)}'
        factors_in_isoform = merged_isoforms[iso_name]['factors']
        EL_sum = list(EL[factors_in_isoform].sum(axis=1))
        EL_merged[iso_name] = EL_sum
    
    return EL_merged

def process_gene(gene, samples, gtf_fh, eta_fh, pvals_fh, H_fh, kw_pvals_fh):
    rds = load_ebpmf_gene(gene)
    if 'strand' not in rds['rds'].keys():
        return None
    strand = rds['rds']['strand'][0]
    if strand == 'plus':
        strand = '+'
    else:
        strand = '-'
        
    junc_file = '/project2/mstephens/cfbuenabadn/gtex-stm/code/junctions.tab.gz'
    coords = rds['rds']['ebpmf_10']['coords']
    junctions_bed = get_junctions_bed(coords, junc_file, gene)
    coords = list(rds['rds']['ebpmf_10']['coords'])
    isoforms = get_isoforms(rds['rds']['ebpmf_10']['train_fit']['EF_smooth'], junctions_bed, coords, binary_fraction=0.5)
    merged_isoforms = merge_isoforms(isoforms)

    gene_gtf = write_gtf(merged_isoforms, gene, strand)
    gtf_fh.write(gene_gtf)

    EL_isoforms = merge_EL(rds['rds']['ebpmf_10']['train_fit']['EL'], merged_isoforms)
    EL_isoforms.index = list(rds['rds']['ebpmf_10']['train_samples'])

        
    tissues = sorted(['Brain_Anterior_cingulate_cortex_BA24',
                      'Brain_Cortex',
                      'Brain_Frontal_Cortex_BA9',
                      'Brain_Putamen_basal_ganglia',
                      'Skin_Not_Sun_Exposed_Suprapubic',
                      'Liver', 'Lung', 'Heart_Atrial_Appendage', 'Muscle_Skeletal', 'Whole_Blood'])
    
    sample_list = []
    for tissue in tissues:
        tissue_samples = list(samples.loc[samples.tissue_id == tissue].index.intersection(EL_isoforms.index))
        sample_list.extend(tissue_samples)

    K = len(merged_isoforms)
    EL_manova = EL_isoforms.loc[sample_list]
    EL_manova['tissue_id'] = list(samples.loc[sample_list, 'tissue_id'])


    eta_row = [gene]
    pvals_row = [gene]

    H_row = [gene]
    kw_pvals_row = [gene]
    if K < 2:
        eta_row.append('NA')
        pvals_row.append('NA')
        H_row.append('NA')
        kw_pvals_row.append('NA')
    else:
    
        for isoform in merged_isoforms.keys():

            try:
                eta_squared, pvals, H, pvals_kw = run_anova(EL_manova, isoform)
            except:
                eta_squared = 'NA'
                pvals = 'NA'
                H = 'NA'
                pvals_kw = 'NA'
            eta_row.append(str(eta_squared))
            pvals_row.append(str(pvals))
            H_row.append(str(H))
            kw_pvals_row.append(str(pvals_kw))

    eta_row = '\t'.join(eta_row) + '\n'
    pvals_row = '\t'.join(pvals_row) + '\n'
    H_row = '\t'.join(H_row) + '\n'
    kw_pvals_row = '\t'.join(kw_pvals_row) + '\n'

    eta_fh.write(eta_row)
    pvals_fh.write(pvals_row)
    H_fh.write(H_row)
    kw_pvals_fh.write(kw_pvals_row)

if __name__ == '__main__':
    with open('ebpmf_models/filtered/tables/snmf.gtf', 'w') as gtf_fh:
        with open('ebpmf_models/filtered/tables/snmf_eta.tab', 'w') as eta_fh:
            with open('ebpmf_models/filtered/tables/snmf_pvals.tab', 'w') as pvals_fh:
                with open('ebpmf_models/filtered/tables/snmf_H.tab', 'w') as H_fh:
                    with open('ebpmf_models/filtered/tables/snmf_kw_pvals.tab', 'w') as kw_pvals_fh:
                        for gene in tqdm(genes, position=0, leave=True):
                            process_gene(gene, samples, gtf_fh, eta_fh, pvals_fh, H_fh, kw_pvals_fh)


    
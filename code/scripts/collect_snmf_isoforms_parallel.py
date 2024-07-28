print('loading modules')
import numpy as np
import pandas as pd
#import rpy2
import rpy2.robjects as ro
#from rpy2.robjects import pandas2ri
#import rpy2.robjects.packages as rpackages
#ebpmf = rpackages.importr('ebpmf.alpha')
import os

print('some more modules')

from matplotlib import pyplot as plt
import seaborn as sns
import tabix
from pybedtools import BedTool
import sys

print('custom-made modules')

sys.path.append('/project2/mstephens/cfbuenabadn/gtex-stm/code/scripts')
from get_isoforms import *
from prepare_counts import *
# from sNMF_plots import *
# from tqdm import tqdm
# import pingouin as pg

print('finished loading modules')

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
        factors_list = merged_isoforms[isoform]['factors']
        factors_description = ':'.join(factors_list)
        iso_start = df.start.astype(int).min()
        iso_end = df.end.astype(int).max()

        gene_start = np.nanmin([gene_start, iso_start])
        gene_end = np.nanmax([gene_end, iso_end])

        iso_start = str(iso_start)
        iso_end = str(iso_end)

        isoform_description = gene_description + f'; transcript_id "{isoform_name}"; factors_id "{factors_description}"'


        isoform_row = [chrom, 's-NMF', 'transcript', iso_start, iso_end, '.', strand, '.', isoform_description]
        isoform_row = '\t'.join(isoform_row)

        gtf_lines.append(isoform_row)

        exon_count = 1
        for idx, row in df.iterrows():
            exon_start = str(int(row.start)) #str(row.start.astype(int))
            exon_end = str(int(row.end)) #str(row.end.astype(int))
            exon_name = f'exon_{str(exon_count)}'
            exon_description = isoform_description + f'; exon_id "{exon_name}"'
            exon_row = [chrom, 's-NMF', 'exon', exon_start, exon_end, '.', strand, '.', exon_description]
            exon_row = '\t'.join(exon_row)
            gtf_lines.append(exon_row)
            exon_count += 1

    gene_start = str(gene_start)
    gene_end = str(gene_end)

    gene_line = [chrom, 's-NMF', 'gene', str(int(float(gene_start))), str(int(float(gene_end))), '.', strand, '.', gene_description]
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
    except:
        print(res)
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

def process_gene(gene, samples, gtf_fh, correct_bias, K):#, eta_fh, pvals_fh, H_fh, kw_pvals_fh):
    rds = load_ebpmf_gene(gene)
    if 'strand' not in rds['rds'].keys():
        return None
    strand = rds['rds']['strand'][0]
    if strand == 'plus':
        strand = '+'
    else:
        strand = '-'
        
    junc_file = '/project2/mstephens/cfbuenabadn/gtex-stm/code/junctions.tab.gz'
    coords = rds['rds'][f'ebpmf_{str(K)}']['coords']
    junctions_bed = get_junctions_bed(coords, junc_file, gene)
    coords_all = list(rds['rds']['coords'])
    coords = list(rds['rds'][f'ebpmf_{str(K)}']['coords'])
    isoforms = get_isoforms(rds['rds'][f'ebpmf_{str(K)}']['train_fit']['EF_smooth'], junctions_bed, coords, coords_all, binary_fraction=0.5, 
                            correct_bias = correct_bias, print_gene = gene)
    merged_isoforms = merge_isoforms(isoforms)

    gene_gtf = write_gtf(merged_isoforms, gene, strand)
    gtf_fh.write(gene_gtf)


if __name__ == '__main__':

    arguments = sys.argv
    
    i = int(arguments[1])
    n = int(arguments[2])
    correction_tag = str(arguments[3])
    K = int(arguments[4])
    
    if correction_tag == '3prime_corrected':
        run_lm = True
    elif correction_tag == 'uncorrected':
        run_lm = False
    else:
        raise Exception('Argument for correcting 3 prime bias should be "3prime_corrected" or "uncorrected"')

    print('Reading list of genes.')
    genes = os.listdir('ebpmf_models/filtered/RDS/')
    genes = [x.split('.')[0] for x in genes]
    print('total genes:')
    print(len(genes))
    print('')
    samples = pd.read_csv('config/samples.tsv', sep='\t', index_col=0)
    print('loaded samples files')

    total_genes = len(genes)

    genes_per_slice = int(total_genes/n)

    slice_start = (i-1)*genes_per_slice
    if i == n:
        slice_end = total_genes + 1
    else:
        slice_end = i*genes_per_slice

    selected_genes = genes[slice_start:slice_end]

    with open(f'ebpmf_models/filtered/snmf_{str(K)}/tables/tmp/snmf.{correction_tag}_{str(i)}.gtf', 'w') as gtf_fh:
        counter = 1
        print('start processing gene')
        selected_genes = list(selected_genes)
        for gene in selected_genes:
            if gene in genes:
                # print(gene)
                if counter % 100 == 0:
                    print(counter)
                counter += 1
                try:
                    process_gene(gene, samples, gtf_fh, correct_bias = run_lm, K=K)#, eta_fh, pvals_fh, H_fh, kw_pvals_fh)
                except:
                    continue
                                    
    print('done!')


    

import numpy as np
import pandas as pd
import gtfparse
import gzip
from pybedtools import BedTool
import tabix

import sys


def get_juncs(bed):
    juncs = pd.DataFrame()
    chrom_cols = list(bed.chrom)[:-1]
    start_cols = list(bed.end.astype(int)+1)[:-1]
    end_cols = list(bed.start.astype(int)-1)[1:]
    juncs['chrom'] = chrom_cols
    juncs['start'] = start_cols
    juncs['end'] = end_cols
    return juncs

def same_intron_chain(ref, query):
    if (len(ref) != len(query)) or (len(query) == 1):
        return False
    else:
        ref_juncs = get_juncs(ref)
        query_juncs = get_juncs(query)
        match_test = np.array(ref_juncs == query_juncs)
        is_chain = np.all(match_test)
        # ref_coords = ref[['start', 'end']].astype(int)
        # query_coords = query[['start', 'end']].astype(int)

        # print((ref_coords - query_coords).start.abs())
        # print((ref_coords - query_coords).end.abs()[:-1].sum())

        # start_diff = (ref_coords - query_coords).start.abs()[1:].sum() == 0
        # end_diff = (ref_coords - query_coords).end.abs()[:-1].sum() == 0

        if is_chain:
            return True
        else:
            return False

def alt_utr(ref, query, strand, prime = 'five'):
    if strand == '-':
        ref = ref[::-1]
        query = query[::-1]

    if prime == 'five':
        i = 0
    elif prime == 'three':
        i = -1
    else:
        raise Exception('invalid prime end')
        
    ref_start = int(ref.iloc[i].start)
    ref_end = int(ref.iloc[i].end)
    ref_exon_len = ref_end - ref_start

    query_start = int(query.iloc[i].start)
    query_end = int(query.iloc[i].end)
    query_exon_len = query_end - query_start

    ratio_ref = ref_exon_len/query_exon_len
    ratio_query = query_exon_len/ref_exon_len
    

    if ((prime == 'five') and (strand == '+')) or ((prime == 'three') and (strand == '-')):
        if ref_end != query_end:
            utr = 'alt.exon'
        elif (np.abs(ref_start - query_start) < 200) and (ratio_ref > 0.75) and (ratio_query > 0.75):
            utr = 'same'
        else:
            utr = 'alt.utr'
    else:
        if ref_start != query_start:
            utr = 'alt.exon'
        elif (np.abs(ref_end - query_end) < 200) and (ratio_ref > 0.75) and (ratio_query > 0.75):
            utr = 'same'
        else:
            utr = 'alt.utr'
    return utr

    
def has_cassette_exon(query, ref):
    query_bed = BedTool.from_dataframe(query)
    ref_juncs = BedTool.from_dataframe(get_juncs(ref))

    cassette_exon = ref_juncs.intersect(query_bed, F=1)

    if len(cassette_exon) >= 1:
        return True
    else:
        return False

def has_missed_exon(query, ref):
    exon_is_missing = has_cassette_exon(ref, query)
    return exon_is_missing

def has_retained_intron(query, ref):
    query_bed = BedTool.from_dataframe(query)
    ref_juncs = BedTool.from_dataframe(get_juncs(ref))

    retained_intron = query_bed.intersect(ref_juncs, F=1)
    if len(retained_intron) >= 1:
        return True
    else:
        return False

def has_novel_intron(query, ref):
    novel_intron =  has_retained_intron(ref, query)
    return novel_intron
                
            
def is_intron_subset(query, ref):
    if (len(query) == 1) or (len(ref)==1):
        return False
    
    query_juncs = BedTool.from_dataframe(get_juncs(query))
    ref_juncs_bed = get_juncs(ref)
    ref_juncs = BedTool.from_dataframe(ref_juncs_bed)

    junc_intersection = query_juncs.intersect(ref_juncs, F=1, f=1).to_dataframe()

    all_juncs_contained = ((len(junc_intersection) == len(query_juncs)) and (len(ref_juncs_bed) > len(junc_intersection)))

    return all_juncs_contained

def is_intron_supraset(query, ref):
    is_supraset = is_intron_subset(ref, query)
    return is_supraset

def is_transcript_subset(query, ref, strand):

    is_chain = same_intron_chain(ref, query)
    is_subset = is_intron_subset(query, ref)

    cassette_exon = has_cassette_exon(query, ref)
    retained_intron = has_retained_intron(query, ref)
    novel_intron = has_novel_intron(query, ref)
    skipped_exon = has_missed_exon(query, ref)

    if (not is_chain) and (not cassette_exon) and (not retained_intron) and (not novel_intron) and (not skipped_exon) and is_subset:
        query_bed = BedTool.from_dataframe(query)
        ref_bed = BedTool.from_dataframe(ref)
        ref_subset = ref_bed.intersect(query_bed, wa=True).to_dataframe()

        alt_5 = alt_utr(ref_subset, query, strand, prime = 'five')
        alt_3 = alt_utr(ref_subset, query, strand, prime = 'three')

        if (alt_5 != 'alt.exon') and (alt_3 != 'alt.exon'):
            return True
            
    return False


def is_transcript_supraset(query, ref, strand):
    is_supraset = is_transcript_subset(ref, query, strand)
    return is_supraset

def is_same_transcript(query, ref, strand):
    if len(query) != len(ref):
        return False
    if len(query) == 1:
        is_chain = True
    else:
        is_chain = same_intron_chain(ref, query)
    alt_5 = alt_utr(ref, query, strand, prime = 'five')
    alt_3 = alt_utr(ref, query, strand, prime = 'three')

    if is_chain and (alt_5 == 'same') and (alt_3 == 'same'):
        return True
    else:
        return False

def annotate_transcript(query, ref, strand):
    is_same = is_same_transcript(query, ref, strand)
    if is_same:
        return 'ref.transcript'
    alt_5 = alt_utr(query, ref, strand, 'five')
    alt_3 = alt_utr(query, ref, strand, 'three')
    out = []
    if alt_5 == 'alt.utr':
        out.append('utr5')
    elif (alt_5 == 'alt.exon') and (len(query) > 1):
        out.append('alt.5exon')
    if alt_3 == 'alt.utr':
        out.append('utr3')
    elif (alt_3 == 'alt.exon') and (len(query) > 1):
        out.append('alt.3exon')
    is_chain = same_intron_chain(query, ref)
    if is_chain:

        out = ['intron.chain'] + out
        
        assert len(out) > 0
        out = '/'.join(out)
        return out
    else:
        cassette_exon = has_cassette_exon(query, ref)
        retained_intron = has_retained_intron(query, ref)
        novel_intron = has_novel_intron(query, ref)
        skipped_exon = has_missed_exon(query, ref)

        if cassette_exon:
            out.append('cassette.exon')
        if retained_intron:
            out.append('retained.intron')
        if novel_intron:
            out.append('novel.intron')
        if skipped_exon:
            out.append('skipped.exon')

        alt_ss_count = count_alt_ss(query, ref)
        if alt_ss_count > 0:
            alt_ss_annot = f'alt.ss.{str(alt_ss_count)}'
            out.append(alt_ss_annot)
    # if len(out) > 0:
        # out = '/'.join(out)
        # return out
    # else:
    if is_transcript_subset(query, ref, strand):
        out = ['transcript.subset'] + out
    elif is_transcript_supraset(query, ref, strand):
        out = ['transcript.supraset'] + out
    if len(out) >= 1:
        out = '/'.join(out)
        return out
    else:
        return ('other')
    

def count_alt_ss(query, ref):

    alt_ss = 0
    
    if (len(query) < 3) or (len(ref) < 3):
        return alt_ss
        
    query_inner_exons = query[1:-1]
    ref_inner_exons = ref[1:-1]

    for idx, row in query_inner_exons.iterrows():
        for idx2, row2 in ref_inner_exons.iterrows():
            if (int(row.start) == int(row2.start)) and (int(row.end) != int(row2.end)):
                alt_ss += 1
            elif (int(row.start) != int(row2.start)) and (int(row.end) == int(row2.end)):
                alt_ss += 1
    return alt_ss
        

def annotate_transcripts_with_ref(query_bed, ref_bed):
    
    strand = list(query_bed.strand)[0]
    assert strand == list(ref_bed.strand)[0]

    bed_cols = ['chrom', 'start', 'end']
    
    ref_transcript_id = get_best_reference(ref_bed)
    ref_transcript_bed = ref_bed.loc[ref_bed.transcript_id == ref_transcript_id, bed_cols].reset_index(drop=True)

    query_transcripts = sorted(query_bed.transcript_id.unique())
    gencode_transcripts = [x for x in sorted(ref_bed.transcript_id.unique()) if x != ref_transcript_id]
    appris_annotation = []
    gencode_annotation = []

    reference_transcript_len = len(ref_transcript_bed)
    query_transcript_len = []

    for query_transcript in query_transcripts:
        query_transcript_bed = query_bed.loc[query_bed.transcript_id == query_transcript, bed_cols].reset_index(drop=True)
        query_transcript_len.append(len(query_transcript_bed))
        query_appris_annotation = annotate_transcript(query_transcript_bed, ref_transcript_bed, strand)
        appris_annotation.append(query_appris_annotation)

        query_gencode_annotation = '.'
        for gencode_transcript in gencode_transcripts:
            gencode_transcript_bed = ref_bed.loc[ref_bed.transcript_id == gencode_transcript, bed_cols].reset_index(drop=True)
            query_gencode = annotate_transcript(query_transcript_bed, gencode_transcript_bed, strand)
            if query_gencode == 'ref.transcript':
                query_gencode_annotation = 'ref.transcript:' + gencode_transcript
                break
            elif 'intron.chain' not in query_gencode_annotation:
                if 'intron.chain' in query_gencode:
                    query_gencode_annotation = 'intron.chain:' + gencode_transcript
                elif 'transcript.subset' in query_gencode:
                    query_genocde_annotation = 'transcript.subset:' + gencode_transcript

        gencode_annotation.append(query_gencode_annotation)

    #appris_annotation = [x + ':' + ref_transcript_id for x in appris_annotation]

    df_out = pd.DataFrame()
    df_out['transcript'] = query_transcripts
    df_out['appris_transcript'] = [ref_transcript_id]*len(query_transcripts)
    df_out['appris_transcript_length'] = [reference_transcript_len]*len(query_transcripts)
    df_out['query_transcript_len'] = query_transcript_len
    df_out['appris_ref'] = appris_annotation
    df_out['gencode_ref'] = gencode_annotation

    return df_out

def get_best_reference(ref_exons):
    if any(ref_exons.appris == 'appris_principal_1'):
        ref_exons = ref_exons.loc[ref_exons.appris == 'appris_principal_1']
    if any(ref_exons.MANE_Select == 'MANE_Select'):
        ref_exons = ref_exons.loc[ref_exons.MANE_Select == 'MANE_Select']
    if any(ref_exons.Ensembl_canonical == 'Ensembl_canonical'):
        ref_exons = ref_exons.loc[ref_exons.Ensembl_canonical == 'Ensembl_canonical']
    if any(ref_exons.basic == 'basic'):
        ref_exons = ref_exons.loc[ref_exons.basic == 'basic']
    if any(ref_exons.transcript_support_level.astype(int) == 1):
        ref_exons = ref_exons.loc[ref_exons.transcript_support_level.astype(int) == 1]
    selected_transcript = list(ref_exons.transcript_id)[0]
    return selected_transcript

if __name__ == '__main__':
    arguments = sys.argv
    gencode_exons_bed = arguments[1]
    snmf_exons_bed = arguments[2]
    # out_prefix = arguments[3]
    i = int(arguments[3])
    n = int(arguments[4])
    correction_tag = str(arguments[5])
    K = str(int(arguments[6]))

    # gencode_exons_bed = 'gencode.v44.primary_assembly.exons.sorted.bed.gz'
    
    gencode_exons = pd.read_csv(gencode_exons_bed, sep='\t', names = ['chrom', 'start', 'end', 'gene_id', 
                                                           'transcript_id', 'strand', 'exon_id', 'transcript_support_level',
                                                           'basic', 'Ensembl_canonical', 'MANE_Select', 'appris', 'transcript_type'])
    
    # snmf_exons_bed = '/project2/mstephens/cfbuenabadn/gtex-stm/code/ebpmf_models/filtered/snmf_10/tables/snmf.exons.sorted.bed.gz'
    
    snmf_exons = pd.read_csv(snmf_exons_bed, sep='\t', names = ['chrom', 'start', 'end', 'gene_id', 
                                                     'transcript_id', 'strand', 'factors', 'exon_id'])
    genes = list(snmf_exons.gene_id.unique())
    total_genes = len(genes)

    genes_per_slice = int(total_genes/n)

    slice_start = (i-1)*genes_per_slice
    if i == n:
        slice_end = total_genes + 1
    else:
        slice_end = i*genes_per_slice

    selected_genes = genes[slice_start:slice_end]

    out_prefix = f'ebpmf_models/filtered/snmf_{K}/tables/tmp/annotated.snmf.{correction_tag}'
    
    annot_out_tab = out_prefix + f'_{str(i)}.tab.gz'
    
    with gzip.open(annot_out_tab, 'wb') as fh:
        write_first_line = True
        counts = 1
        
        for gene in selected_genes:

            try:
            
                snmf_bed = snmf_exons.loc[snmf_exons.gene_id == gene]
                gencode_bed = gencode_exons.loc[gencode_exons.gene_id == gene]
                annot = annotate_transcripts_with_ref(snmf_bed, gencode_bed)
                if write_first_line:
                    first_line = ('\t'.join(list(annot.columns)) + '\n').encode()
                    fh.write(first_line)
                    write_first_line = False
                for idx, row in annot.iterrows():
                    line = ('\t'.join(list(row.astype(str))) + '\n').encode()
                    fh.write(line)
            except:
                continue
        
            # if counts % 10 == 0:
            #     print(counts)
            # counts += 1
                
            
            

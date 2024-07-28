def get_bed_per_tissue_15(wildcards):
    bed_samples = []
    tissue_15 = tissue_list + new_tissue
    for tissue in tissue_15:
        tissue_samples = gtex_samples.loc[gtex_samples.tissue_id == tissue].index
        file_name = "coverage/bed/" + tissue + "/{IndID}.bed.gz"
        bed_samples.extend(expand(file_name, IndID=tissue_samples))
    return bed_samples
    
def get_tbi_per_tissue_15(wildcards):
    tbi_samples = []
    tissue_15 = tissue_list + new_tissue
    for tissue in tissue_15:
        tissue_samples = gtex_samples.loc[gtex_samples.tissue_id == tissue].index
        file_name = "coverage/bed/" + tissue + "/{IndID}.bed.gz.tbi"
        tbi_samples.extend(expand(file_name, IndID=tissue_samples))
    return tbi_samples

rule MakeGeneCounts_15:
    input:
        get_bed_per_tissue_15,
        get_tbi_per_tissue_15
    output:
        "coverage/counts_total/{gene}.csv.gz"
    resources:
        mem_mb = 24000,
    log:
        "/scratch/midway3/cnajar/logs/counts/tissues.{gene}.log"
    wildcard_constraints:
        gene = '|'.join(selected_genes.gene),
    params:
        tissues = ' '.join(tissue_list)
    shell:
        """
        python scripts/prepare_counts.py {wildcards.gene} {params.tissues} &> {log}
        """

rule FilterCountTables:
    input:
        "coverage/counts_total/{gene}.csv.gz"
        #expand("coverage/counts/{Tissue}/{{gene}}.csv.gz", Tissue = tissue_list)
    output:
        "coverage/counts_filtered/{gene}.csv.gz",
        "coverage/counts_filtered_stats/{gene}.stats"
        #"coverage/count_tables/{gene}.tab.gz"
    log:
        '/scratch/midway3/cnajar/logs/make_counts/{gene}.log'
    resources:
        mem_mb = 68000,
    shell:
        """
        python scripts/filter_counts.py {wildcards.gene} &> {log}
        """


rule MakeGeneCounts_WholeGene:
    input:
        get_bed_per_tissue,
        get_tbi_per_tissue
    output:
        "coverage/counts_whole_gene/{gene}.csv.gz"
    resources:
        mem_mb = 24000,
    log:
        "/scratch/midway3/cnajar/logs/counts_whole_gene/tissues.{gene}.log"
    wildcard_constraints:
        gene = '|'.join(selected_genes.gene),
        #Tissue = '|'.join(tissue_list)
    params:
        tissues = ' '.join(tissue_list)
    shell:
        """
        python scripts/prepare_counts_whole_gene.py {wildcards.gene} {params.tissues} &> {log}
        """


rule run_snmf:
    input:
        "coverage/15_tissues/counts_{Pass}/{gene}.csv.gz"
    output:
        'ebpmf_models/15_tissues/{Pass}/RDS/{gene}.rds',
    log:
        '/scratch/midway3/cnajar/logs/ebpmf_run/{gene}.{Pass}.log'
    resources:
        mem_mb = 48000,
    wildcard_constraints:
        gene = '|'.join(list(selected_genes.gene)),
        Pass = 'whole_gene|filtered'
    params:
        strand = get_strand
    shell:
        """
        {config[Rscript]} scripts/run_snmf.R {wildcards.gene} {params.strand} {wildcards.Pass} &> {log}
        """


rule ebpmf_run_single_tissue:
    input:
        "coverage/15_tissues/counts_filtered/{gene}.csv.gz"
    output:
        'ebpmf_models/single_tissue/RDS/{gene}.rds',
    log:
        'logs/ebpmf_run/{gene}.single_tissue.log'
    resources:
        mem_mb = 24000,
    wildcard_constraints:
        gene = '|'.join(list(selected_genes.gene)),
    params:
        strand = get_strand
    shell:
        """
        ({config[Rscript]} scripts/run_ebpmf_single_tissue.R {wildcards.gene} {params.strand}) &> {log}
        """
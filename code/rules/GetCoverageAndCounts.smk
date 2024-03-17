def much_more_mem_after_first_attempt(wildcards, attempt):
    if int(attempt) == 1:
        return 24000
    else:
        return 62000


#rule GetBigWigs:
#    input:
#        bam = "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{Tissue}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
#        bai = "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-#download/{Tissue}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.bai"
#    output:
#        "coverage/bigwigs/{Tissue}/{IndID}.bw"
#    resources:
#        mem_mb = 42000,
#    log:
#        "logs/bigwigs/{Tissue}.{IndID}.log"
#    shell:
#        """
#        (bamCoverage -b {input.bam} -o {output}) &> {log}
#        """


rule GetBedsAndTabix:
    input:
        bam = "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{Tissue}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        bai = "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{Tissue}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.bai"
    output:
        bed = "coverage/bed/{Tissue}/{IndID}.bed.gz",
        tabix = "coverage/bed/{Tissue}/{IndID}.bed.gz.tbi",
    resources:
        mem_mb = much_more_mem_after_first_attempt
    wildcard_constraints:
        Tissue = "|".join(tissue_list)
    log:
        "logs/bed_and_tabix/{Tissue}.{IndID}.log"
    shell:
        """
        (bedtools genomecov -5 -bga -ibam {input.bam} | bgzip > {output.bed}) &> {log};
        (tabix -p bed {output.bed}) &>> {log}
        """

def get_bed_per_tissue(wildcards):
    bed_samples = []
    for tissue in tissue_list:
        tissue_samples = gtex_samples.loc[gtex_samples.tissue_id == tissue].index
        file_name = "coverage/bed/" + tissue + "/{IndID}.bed.gz"
        bed_samples.extend(expand(file_name, IndID=tissue_samples))
    return bed_samples
    
def get_tbi_per_tissue(wildcards):
    tbi_samples = []
    for tissue in tissue_list:
        tissue_samples = gtex_samples.loc[gtex_samples.tissue_id == tissue].index
        file_name = "coverage/bed/" + tissue + "/{IndID}.bed.gz.tbi"
        tbi_samples.extend(expand(file_name, IndID=tissue_samples))
    return tbi_samples

#tissue_samples = gtex_samples.loc[gtex_samples.tissue_id == wildcards.Tissue].index
#tbi_samples = expand("coverage/bed/{{Tissue}}/{IndID}.bed.gz.tbi", IndID=tissue_samples)
#return tbi_samples
    
def get_gene_coords(wildcards):
    gene_bed = selected_genes.loc[selected_genes.gene==wildcards.gene]
    gene_chrom = gene_bed.chrom
    gene_start = str(int(gene_bed.start) - 50)
    gene_end = str(int(gene_bed.end) + 50)
    coords = ' '.join([gene_chrom, gene_start, gene_end])
    return coords

rule MakeGeneCounts:
    input:
        get_bed_per_tissue,
        get_tbi_per_tissue
    output:
        temp("coverage/counts_total/{gene}.csv.gz")
        #temp(expand("coverage/counts/{Tissue}/{{gene}}.csv.gz", Tissue=tissue_list))
    resources:
        mem_mb = 24000,
    log:
        "/scratch/midway3/cnajar/logs/counts/tissues.{gene}.log"
    wildcard_constraints:
        gene = '|'.join(selected_genes.gene),
        #Tissue = '|'.join(tissue_list)
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
    
    
        

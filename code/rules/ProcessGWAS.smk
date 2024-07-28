rule format_GWAS_summary_statistics:
    output:
        expand('gwas/temp/GRCh37/{trait}.bed', trait = gwas_annot.loc[gwas_annot.genome_assembly=='GRCh37'].trait),
        expand('gwas/temp/GRCh38/{trait}.bed', trait = gwas_annot.loc[gwas_annot.genome_assembly=='GRCh38'].trait)
    log:
        'logs/gwas/formatting_summary_stats.log'
    resources:
        mem_mb = 12000
    shell:
        """
        bash scripts/process_gwas.sh &> {log}
        """

# for I in $(grep -v figshare ../config/gwas.txt | awk -F'\t' '{print $6}'); do wget $I; done

rule liftover_Gwas_stats:
    """
    Convert hg19 bed of summary stats to hg38_summarystats. For convenience with downstream rules, make sure the input bed has Pvalues in column4. Other summary stats can be in later columns
    """
    input:
        bed = "gwas/temp/GRCh37/{trait}.bed",
        chain = "Annotations/hg19ToHg38.over.chain.gz"
    output:
        "gwas/temp/GRCh38/{trait}.bed"
    shadow: "shallow"
    wildcard_constraints:
        trait = '|'.join(list(gwas_annot.loc[gwas_annot.genome_assembly=='GRCh37'].trait))
    shell:
        """
        CrossMap.py bed {input.chain} {input.bed} gwas/temp/GRCh38/{wildcards.trait}.bed
        """


def get_standard_conversion(wildcards):
    summary_statistics = gwas_annot.loc[gwas_annot.trait == wildcards.trait, 'summary_statistics'].iloc[0]
    return summary_statistics

rule StandardizeSummaryStatistics:
    input:
        "gwas/temp/GRCh38/{trait}.bed"
    output:
        "gwas/temp/GRCh38/{trait}.beta_se.bed"
    log:
        'logs/gwas/standardize.summary_stats.{trait}.log'
    resources:
        mem_mb = 12000
    params:
        get_standard_conversion
    wildcard_constraints:
        trait = '|'.join(list(gwas_annot.trait))
    shell:
        """
        Rscript scripts/StandardizeGwasStats.R {input} {params} {output} &> {log}
        """



rule SortCompressAndIndex:
    input:
        "gwas/temp/GRCh38/{trait}.beta_se.bed"
    output:
        bed = temp("gwas/hg38_summary_stats/{trait}.bed"),
        bedgz = "gwas/hg38_summary_stats/{trait}.bed.gz",
        tbi = "gwas/hg38_summary_stats/{trait}.bed.gz.tbi"
    log:
        'logs/gwas/standardize.sort_compress_index.{trait}.log'
    resources:
        mem_mb = 58000
    wildcard_constraints:
        trait = '|'.join(list(gwas_annot.trait))
    shell:
        """
        (head -n1 {input} > {output.bed}) &> {log};
        (tail -n+2 {input} | sort -k 1,1 -k2,2n - >> {output.bed}) &> {log};
        (bgzip /dev/stdin -c {output.bed} > {output.bedgz}) &> {log};
        (tabix -p bed {output.bedgz}) &> {log};
        """

def GetWindowSize(wildcards):
    if wildcards.WindowSize == 'Window1M':
        return 500000
    elif wildcards.WindowSize == 'Window100K':
        return 50000
    elif wildcards.WindowSize == 'Window10K':
        return 5000

rule GetGWAS_LeadSnpWindows:
    """
    output bed of 1MB window surrounding lead genome-wide signficant autosomal
    snps for each gwas. Exclude blacklistregions (ie MHC). This rule only works
    for summary stats in bed format with Pvalue in column4. Exit status 1 if
    output is empty
    """
    input:
        summarystats = "gwas/hg38_summary_stats/{trait}.bed.gz",
        blacklistregions = "../data/MHC.hg38.bed",
        chromsizes = "Annotations/GRCh38.primary_assembly.genome.fa.fai"
    output:
        signif_loci = "gwas/leadSnps/{trait}.bed"
    params:
        PvalThreshold = "5e-8"
    log:
        "logs/GetGWAS_LeadSnpWindows_NonGwasCatalog/{trait}.log"
    wildcard_constraints:
        trait = '|'.join(list(gwas_annot.trait))
    resources:
        mem_mb = 58000
    shell:
        """
        (python scripts/GetGWASLeadVariantWindowsFromBed.py {input.summarystats} /dev/stdout {params.PvalThreshold} | awk -F'\\t' -v OFS='\\t' '$1~/chr[0-9]+/ {{print $1, $2, $2, $1"_"$2"_N_N_{wildcards.trait}" }}' | bedtools slop -i - -g {input.chromsizes} -b 500000 | bedtools sort -i - | bedtools intersect -a - -b {input.blacklistregions} -wa -sorted -v > {output} ) &> {log}
        [[ -s {output.signif_loci} ]]
        """

rule collect_grch38:
    input:
        expand("gwas/leadSnps/{trait}.bed", trait = gwas_annot.trait)



######

rule ConcatGwasLeadSnpWindows:
    input:
        expand(
            "gwas/leadSnps/{trait}.bed", trait=gwas_traits
        ),
    output:
        "gwas/LeadSnpWindows.bed"
    log:
        "logs/ConcatGwasLeadSnpWindows.log"
    shell:
        """
        cat {input} | bedtools sort -i - > {output}
        """


rule GwasBedStatsAtWindows:
    """
    For coloc, gather summary stats in window centered on lead SNPs
    """
    input:
        signif_loci = "gwas/leadSnps/{trait}.bed",
        bed = "gwas/hg38_summary_stats/{trait}.bed.gz",
        tbi = "gwas/hg38_summary_stats/{trait}.bed.gz.tbi",
    log:
        "logs/GwasBedStatsAtWindows/{trait}.log"
    output:
        stats = "gwas/StatsForColoc/{trait}.standardized.txt.gz"
    wildcard_constraints:
        accession = '|'.join([x for x in list(gwas_annot.trait) if x in gwas_traits])
    shell:
        """
        (tabix -h -R {input.signif_loci} {input.bed} | sort -k1,1 -k2,2n | bedtools intersect -sorted -a - -b {input.signif_loci} -wo | awk -F'\\t' -v OFS='\\t'  '{{print $12, $1, $2, $7, $8, $6}}'  | awk -F'[:\\t]' 'BEGIN{{print "loci\\tchrom\\tstart\\tbeta\\tSE\\tA1\\tA2"}} {{print $1, $2, $3, $4, $5, $8, $9}}' OFS='\\t' | gzip - > {output.stats} ) &> {log}
        """

rule CollectLeadSnps:
    input:
        expand("gwas/StatsForColoc/{trait}.standardized.txt.gz", trait=[x for x in list(gwas_annot.trait) if x in gwas_traits]),
        "gwas/LeadSnpWindows.bed"
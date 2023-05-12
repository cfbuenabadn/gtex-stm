rule makeGene_bed:
    input:
        "Annotations/genes.tab.gz",
    output:
        "coverage/tmp/{gene}.bed"
    wildcard_constraints:
        gene = "|".join(genes)
    log:
        "logs/genebed.{gene}_bed.log"
    shell:
        """
        python scripts/makebed.py {wildcards.gene} &> {log}
        """

def much_more_mem_after_first_attempt(wildcards, attempt):
    if int(attempt) == 1:
        return 24000
    else:
        return 72000

rule MakeBedFromBam:
    input:
        bam = "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{Tissue}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        bai = "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{Tissue}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.bai"
    output:
        "coverage/samples/{Tissue}/{IndID}.bed.gz"
    log:
        "logs/bamcoverage.{Tissue}.{IndID}.log"
    resources:
        mem_mb = much_more_mem_after_first_attempt
    wildcard_constraints:
        gene = "|".join(genes),
        Tissue = '|'.join(tissue_list)
    shell:
        """
        (bedtools genomecov -5 -bga -ibam {input.bam} | bedtools sort -i - | gzip - > {output}) &> {log}
        """
        
rule GetGeneBed:
    input:
        bed = "coverage/samples/{Tissue}/{IndID}.bed.gz",
        gene_bed = "coverage/tmp/{Gene}.bed"
    output:
        temp("coverage/bed/{Gene}/{Tissue}.{IndID}.bed.gz")
    #log:
    #    "logs/genecoverage.{Gene}.{Tissue}.{IndID}.log"
    resources:
        mem_mb = much_more_mem_after_first_attempt
    wildcard_constraints:
        Gene = "|".join(genes),
        Tissue = '|'.join(tissue_list)
    shell:
        """
        bedtools intersect -a {input.bed} -b {input.gene_bed} | bedtools sort -i - | gzip - > {output}  #&> {log}
        """
        
def GetCountsBedsPerTissue(wildcards):
    tissue_samples = gtex_samples.loc[gtex_samples.tissue_id == wildcards.Tissue].index
    bed_samples = expand("coverage/bed/{gene}/{Tissue}.{{IndID}}.bed.gz".format(gene=wildcards.gene, Tissue=wildcards.Tissue), IndID = tissue_samples)
    return bed_samples

rule GetCountsGene_generalized:
    input:
        GetCountsBedsPerTissue
    output:
        "Counts/{Tissue}/{gene}.Counts.csv.gz"
    log:
        "logs/{Tissue}.{gene}.Counts.log"
    wildcard_constraints:
        Gene = "|".join(genes),
        Tissue = '|'.join(tissue_list)
    shell:
        """
        python scripts/getCountsTable2.py --output {output} {input} &> {log}
        """
        

    
    
    
   

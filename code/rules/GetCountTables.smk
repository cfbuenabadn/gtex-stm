rule GetBamCoverage:
    input:
        bam = "/project2/yangili1/cfbuenabadn/gtex-stm/code/gtex-download/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        bai = "/project2/yangili1/cfbuenabadn/gtex-stm/code/gtex-download/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.bai"
    output:
        gCoverage = "BamCoverage/{IndID}.coverage.bed"
    log:
        "logs/coverage.{IndID}.log"
    shell:
        """
        bedtools genomecov -5 -bga -ibam {input.bam} > {output.gCoverage}
        """
            
rule IntersectCoverageWithGene:
    input:
        gCoverage = "BamCoverage/{IndID}.coverage.bed",
        gene = "Annotations/genes.bed"
    output:
        "BamTranscriptome/{IndID}.transcriptome.bed"
    log:
        "logs/{IndID}.transcriptome.log"
    shell:
        """
        bedtools intersect -wao -a {input.gene} -b {input.gCoverage} > {output}
        """
        
rule GetCoverageGene:
    input:
        "BamTranscriptome/{IndID}.transcriptome.bed"
    output:
        "CoverageGene/{gene}/{IndID}.coverage.bed"
    log:
        "logs/{gene}.{IndID}.coverage.log"
    shell:
        """
        awk -F'\t' '$4=="{wildcards.gene}"' {input} > {output}
        """
        
rule SortBedFiles:
    input:
        "CoverageGene/{gene}/{IndID}.coverage.bed"
    output:
        "CoverageGene/{gene}/{IndID}.coverage.sorted.bed"
    log:
        "logs/{gene}.{IndID}.coverage.sort.log"
    shell:
        """
        sort -k 6 -n {input} > {output}
        """
        
rule GetCountsGene:
    input:
       lambda wildcards: expand("CoverageGene/{{gene}}/{IndID}.coverage.sorted.bed", IndID = gtex_samples)
    output:
        "Counts/{gene}.Counts.csv.gz"
    log:
        "logs/{gene}.Counts.log"
    shell:
        """
        python scripts/getCountsTable.py --output {output} --gene_dir CoverageGene/{wildcards.gene}/ &> {log}
        """
        
        
            
            
            

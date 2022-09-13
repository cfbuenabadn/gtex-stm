rule GetBamCoverage:
    input:
        bam = "/project2/yangili1/cfbuenabadn/gtex-stm/code/gtex-download/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        bai = "/project2/yangili1/cfbuenabadn/gtex-stm/code/gtex-download/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.bai"
    output:
        gCoverage = temp("BamCoverage/{IndID}.coverage.bed")
    log:
        "logs/coverage.{IndID}.log"
    shell:
        """
        bedtools genomecov -5 -bga -ibam {input.bam} > {output.gCoverage}
        """
            
rule IntersectCoverageWithGene:
    input:
        gCoverage = "BamCoverage/{IndID}.coverage.bed",
        genes = "Annotations/genes.bed"
    output:
        "BamTranscriptome/{IndID}.transcriptome.bed.gz"
    log:
        "logs/{IndID}.transcriptome.log"
    resources:
        mem_mb = 58000
    shell:
        """
        bedtools intersect -wao -a {input.genes} -b {input.gCoverage} > BamTranscriptome/{wildcards.IndID}.transcriptome.bed;
        gzip BamTranscriptome/{wildcards.IndID}.transcriptome.bed
        """
        
rule GetCoverageGene:
    input:
        lambda wildcards: expand("BamTranscriptome/{IndID}.transcriptome.bed.gz", IndID = gtex_samples)
    output:
        ["CoverageGene/{gene}/" + IndID + ".coverage.sorted.bed.gz" for IndID in gtex_samples]
        #f"CoverageGene/{{gene}}/{gtex_samples}.coverage.sorted.bed.gz"
        #"CoverageGene/{gene}/token.txt"
        #lambda wildcards: expand("CoverageGene/{{gene}}/{IndID}.coverage.sorted.bed.gz", IndID = gtex_samples)
    log:
        "logs/{gene}.coverage.log"
    resources:
        mem_mb = 58000
    shell:
        """
        for IndID in $(ls BamTranscriptome/ | awk -F'.' '{{print $1}}')
        do
          zcat BamTranscriptome/${{IndID}}.transcriptome.bed.gz | awk -F'\\t' '$4=="{wildcards.gene}"' - > CoverageGene/{wildcards.gene}/${{IndID}}.coverage.bed
          sort -k 6 -n CoverageGene/{wildcards.gene}/${{IndID}}.coverage.bed > CoverageGene/{wildcards.gene}/${{IndID}}.coverage.sorted.bed
          gzip CoverageGene/{wildcards.gene}/${{IndID}}.coverage.sorted.bed
          rm CoverageGene/{wildcards.gene}/${{IndID}}.coverage.bed
        done
        """


rule GetCountsGene:
    input:
        lambda wildcards: expand("CoverageGene/{{gene}}/{IndID}.coverage.sorted.bed.gz", IndID = gtex_samples)
    output:
        "Counts/{gene}.Counts.csv.gz"
    log:
        "logs/{gene}.Counts.log"
    wildcard_constraints:
        gene = '|'.join(genes)
    shell:
        """
        python scripts/getCountsTable.py --output {output} {input} &> {log}
        """
            
            
rule GetCountsGene_NoTissue:
    input:
        "Counts/{gene}.Counts.csv.gz",
        "../data/sample.tsv",
        "../data/participant.tsv"
    output:
        "Counts/{gene}.Counts.no_{tissue}.csv.gz",
        "Counts/{gene}.Counts.one_of_{tissue}.csv.gz",
        "Counts/{gene}.Counts.minus_one_of_{tissue}.csv.gz"
    wildcard_constraints:
        gene = '|'.join(genes)
    log:
        "logs/{gene}.Counts_remove_{tissue}.log"
    shell:
        """
        python scripts/getCountsTableNoTissue.py --gene {wildcards.gene} --tissue {wildcards.tissue} &> {log}
        """
      
      
       

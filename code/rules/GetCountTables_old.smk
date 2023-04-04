#bamCoverage -b SMTN.bam -o cov.bg -of bedgraph -bs 1 -r chr22:31064055:31113807
#bedtools genomecov -5 -bga -ibam SMTN.bam > SMTN.bed

import pandas as pd

def GetGeneCoords(wildcards):
    df = pd.read_csv('../data/genes.tab.gz', names = ['ensembl_id', 'coord'], index_col=0, sep='\t')
    gene_coord = df.loc[wildcards.gene].coord
    return gene_coord

rule MakeGeneBam:
    input:
        gene_tab = "Annotations/genes.tab.gz",
        bam = "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{Tissue}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        bai = "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{Tissue}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.bai"
    output:
        bam = temp("coverage/BAMs/{gene}/{Tissue}.{IndID}.bam"),
        bai = temp("coverage/BAMs/{gene}/{Tissue}.{IndID}.bam.bai")
    wildcard_constraints:
        gene = "|".join(genes)
    log:
        "logs/makegenebam.{gene}.{Tissue}.{IndID}.log"
    params:
        GetGeneCoords
    shell:
        """
        (samtools view -b {input.bam} "{params}" > {output.bam}) &> {log};
        (samtools index {output.bam}) &> {log}
        """

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
        return 64000

rule GetGeneCoverage:
    input:
        bam = "coverage/BAMs/{gene}/{Tissue}.{IndID}.bam",
        bai = "coverage/BAMs/{gene}/{Tissue}.{IndID}.bam.bai",
        bed = "coverage/tmp/{gene}.bed"
    output:
        temp("coverage/bed/{gene}/{Tissue}.{IndID}.bed.gz"),
    log:
        "logs/genecoverage.{gene}.{Tissue}.{IndID}.log"
    resources:
        mem_mb = much_more_mem_after_first_attempt
    wildcard_constraints:
        gene = "|".join(genes),
        Tissue = '|'.join(tissue_list) #['Brain_Cortex', 'Muscle_Skeletal', 'Liver', 'Whole_Blood', 'Lung'])
    shell:
        """
        (bedtools genomecov -5 -bga -ibam {input.bam} | bedtools intersect -a - -b {input.bed} | bedtools sort -i - | gzip - > {output}) &> {log}
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
        gene = "|".join(genes),
        Tissue = '|'.join(tissue_list) #['Brain_Cortex', 'Muscle_Skeletal', 'Liver', 'Whole_Blood', 'Lung'])
    shell:
        """
        python scripts/getCountsTable2.py --output {output} {input} &> {log}
        """



#rule GetCountsGene_BrainCortex2:
#    input:
#        lambda wildcards: expand("coverage/bed/{{gene}}/Brain_Cortex.{IndID}.bed.gz", IndID = brain_cortex_samples)
#    output:
#        "Counts/Muscle_Skeletal/{gene}.Counts2.csv.gz"
#    log:
#        "logs/Muscle_Skeletal/{gene}.Counts.log"
#    wildcard_constraints:
#        gene = 'SMTN'
#    shell:
#        """
#        python scripts/getCountsTable2.py --output {output} {input} &> {log}
#        """



#######################

#rule GetBamCoverage:
#    input:
#        bam = "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{Tissue}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
#        vai = "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{Tissue}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.bai",
#        #bam = "/project2/yangili1/cfbuenabadn/gtex-stm/code/gtex-download/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
#        #bai = "/project2/yangili1/cfbuenabadn/gtex-stm/code/gtex-download/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.bai"
#    output:
#        gCoverage = temp("BamCoverage/{Tissue}/{IndID}.coverage.bed")
#    log:
#        "logs/coverage.{Tissue}.{IndID}.log"
#    shell:
#        """
#        (bedtools genomecov -5 -bga -ibam {input.bam} > {output.gCoverage}) &> {log}
#        """
            
#rule IntersectCoverageWithGene:
#    input:
#        gCoverage = "BamCoverage/{Tissue}/{IndID}.coverage.bed",
#        genes = "Annotations/genes.bed"
#    output:
#        "BamTranscriptome/{Tissue}/{IndID}.transcriptome.bed.gz"
#    log:
#        "logs/{Tissue}.{IndID}.transcriptome.log"
#    resources:
#        mem_mb = 12000
#    shell:
#        """
#        (bedtools intersect -wao -a {input.genes} -b {input.gCoverage} | gzip - > {output}) &> {log};
#        #gzip BamTranscriptome/{wildcards.IndID}.transcriptome.bed
#        """
        
#rule GetCoverageGene_BrainCortex:
#    input:
#        lambda wildcards: expand("BamTranscriptome/Brain_Cortex/{IndID}.transcriptome.bed.gz", IndID = brain_cortex_samples)
#    output:
#        ["CoverageGene/{gene}/" + IndID + ".coverage.sorted.bed.gz" for IndID in brain_cortex_samples]
#    params:
#        tissue = 'Brain_Cortex',
#        samples = ' '.join(brain_cortex_samples)
#    log:
#        "logs/{gene}.coverage.log"
#    resources:
#        mem_mb = 58000
#    shell:
#        """
#        for IndID in {params.samples}
#        do
#          (zcat BamTranscriptome/{params.tissue}/${{IndID}}.transcriptome.bed.gz | awk -F'\\t' '$4=="{wildcards.gene}"' - | sort -k 6 -n - | gzip - > CoverageGene/{wildcards.gene}/${{IndID}}.coverage.sorted.bed.gz) &> {log}
#        done
#        """
#
#use rule GetCoverageGene_BrainCortex as GetCoverageGene_MuscleSkeletal with:
#    input:
#        lambda wildcards: expand("BamTranscriptome/Muscle_Skeletal/{IndID}.transcriptome.bed.gz", IndID = muscle_skeletal_samples)
#    output:
#        ["CoverageGene/{gene}/" + IndID + ".coverage.sorted.bed.gz" for IndID in muscle_skeletal_samples]
#    params:
#        tissue = 'Muscle_Skeletal',
#        samples = ' '.join(muscle_skeletal_samples)

#rule GetCountsGene_BrainCortex:
#    input:
#        lambda wildcards: expand("CoverageGene/Brain_Cortex/{{gene}}/{IndID}.coverage.sorted.bed.gz", IndID = brain_cortex_samples)
#    output:
#        "Counts/Brain_Cortex/{gene}.Counts.csv.gz"
#    log:
#        "logs/Brain_Cortex/{gene}.Counts.log"
#    wildcard_constraints:
#        gene = '|'.join(genes)
#    shell:
#        """
#        python scripts/getCountsTable.py --output {output} {input} &> {log}
#        """

#use rule GetCountsGene_BrainCortex as GetCountsGene_MuscleSkeletal with:
#    input:
#        lambda wildcards: expand("CoverageGene/Muscle_Skeletal/{{gene}}/{IndID}.coverage.sorted.bed.gz", IndID = muscle_skeletal_samples)
#    output:
#        "Counts/Muscle_Skeletal/{gene}.Counts.csv.gz"
#    log:
#        "logs/Muscle_Skeletal/{gene}.Counts.log"
            
            
#rule GetCountsGene_NoTissue:
#    input:
#        "Counts/{gene}.Counts.csv.gz",
#        "../data/sample.tsv",
#        "../data/participant.tsv"
#    output:
#        "Counts/{gene}.Counts.no_{tissue}.csv.gz",
#        "Counts/{gene}.Counts.one_of_{tissue}.csv.gz",
#        "Counts/{gene}.Counts.minus_one_of_{tissue}.csv.gz"
#    wildcard_constraints:
#        gene = '|'.join(genes)
#    log:
#        "logs/{gene}.Counts_remove_{tissue}.log"
#    shell:
#        """
#        python scripts/getCountsTableNoTissue.py --gene {wildcards.gene} --tissue {wildcards.tissue} &> {log}
#        """
      
      
       

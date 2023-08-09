def much_more_mem_after_first_attempt(wildcards, attempt):
    if int(attempt) == 1:
        return 24000
    else:
        return 62000


rule GetBigWigs:
    input:
        bam = "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{Tissue}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        bai = "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{Tissue}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.bai"
    output:
        "coverage/bigwigs/{Tissue}/{IndID}.bw"
    resources:
        mem_mb = 42000,
    log:
        "logs/bigwigs/{Tissue}.{IndID}.log"
    shell:
        """
        (bamCoverage -b {input.bam} -o {output}) &> {log}
        """


rule GetBedsAndTabix:
    input:
        bam = "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{Tissue}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        bai = "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{Tissue}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.bai"
    output:
        bed = "coverage/bed/{Tissue}/{IndID}.bed.bgz",
        tabix = "coverage/bed/{Tissue}/{IndID}.bed.bgz.tbi",
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


rule featureCounts:
    input:
        bam = GetTissueBAM,
        bai = GetTissueBai,
        annotations = "Annotations/gencode.v34.primary_assembly.annotation.gtf"
    output:
        "featureCounts/{Tissue}/Counts.txt"
    threads:
        2
    wildcard_constraints:
        Tissue = "|".join(tissue_list)
    resources:
        mem = 42000,
        cpus_per_node = 2,
    log:
        "logs/featureCounts/{Tissue}.log"
    shell:
        """
        featureCounts -p -T {threads} --ignoreDup --primary -a {input.annotations} -o {output} {input.bam} &> {log}
        """


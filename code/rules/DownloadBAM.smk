rule DownloadFromGTEX:
    input:
        "../data/file-manifest.json",
        "../data/gen3-client"
    output:
        expand(
        "/project2/yangili1/cfbuenabadn/gtex-stm/code/gtex-download/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam", 
        IndID = gtex_samples
        )
    shell:
        """
        ../data/gen3-client download-multiple --profile=AnVIL --manifest=../data/file-manifest.json --download-path=/project2/yangili1/cfbuenabadn/gtex-stm/code/gtex-download/bams/ --protocol=s3
        """

rule GetIndex:
    input:
        bam = "/project2/yangili1/cfbuenabadn/gtex-stm/code/gtex-download/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam"
    output:
        bai = "/project2/yangili1/cfbuenabadn/gtex-stm/code/gtex-download/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.bai"
    log:
        "logs/bam_idx.{IndID}.log"
    shell:
        """
        samtools index {input.bam} {output.bai} &> {log}
        """

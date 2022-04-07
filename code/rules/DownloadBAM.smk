rule DownloadFromGTEX:
    input:
        #manifest = "../data/file-manifest.json",
        #manifest = "../data/file-manifest_20.json",
        manifest = "../data/file-manifest_50.json",
        #manifest = "../data/file-manifest_20-balanced.json",
        #manifest = "../data/file-manifest_30-1.json",
        client = "../data/gen3-client"
    output:
        temp(
        expand(
        "/project2/yangili1/cfbuenabadn/gtex-stm/code/gtex-download/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam", 
        IndID = gtex_samples
        #IndID = gtex_samples_50
        #IndID = gtex_samples_balanced
        #IndID = gtex_samples_30_1
        )
        )
    shell:
        """
        {input.client} download-multiple --profile=AnVIL --manifest={input.manifest} --download-path=/project2/yangili1/cfbuenabadn/gtex-stm/code/gtex-download/bams/ --protocol=s3
        """


rule GetIndex:
    input:
        bam = "/project2/yangili1/cfbuenabadn/gtex-stm/code/gtex-download/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam"
    output:
        bai = temp("/project2/yangili1/cfbuenabadn/gtex-stm/code/gtex-download/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.bai")
    log:
        "logs/bam_idx.{IndID}.log"
    shell:
        """
        samtools index {input.bam} {output.bai} &> {log}
        """

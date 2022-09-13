#rule DownloadFromGTEX:
#    input:
#        manifest = "../data/file-manifest_50.json",
#        client = "../data/gen3-client"
#    output:
#        temp(
#        expand(
#        "/project2/yangili1/cfbuenabadn/gtex-stm/code/gtex-download/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam", 
#        IndID = gtex_samples
#        )
#        )
#    shell:
#        """
#        {input.client} download-multiple --profile=AnVIL --manifest={input.manifest} --download-path=/project2/yangili1/cfbuenabadn/gtex-stm/code/gtex-download/bams/ --protocol=s3
#        """

rule MakeFileManifestJson:
    input:
        "../data/file-manifest.json",
        '../data/sample.tsv',
        '../data/participant.tsv'
    output:
        "gtex-download/file_manifest/TrainTest_2tissues.json",
        "gtex-download/file_manifest/ThreeSamplesFemale.json",
        "gtex-download/file_manifest/ThreeSamplesMale.json",
    log:
        'logs/makefilemanifest.json'
    shell:
        """
        python scripts/get_file_manifest_brain_muscle.py &> {log}
        """

def GetFileManifest(wildcards):
    manifest = "gtex_download/file_manifest/{SampleSet}.json".format(SampleSet = wildcards.SampleSet)
    return manifest

rule DownloadFromGTEX_ThreeSamplesFemale:
    input:
        manifest = "gtex-download/file_manifest/ThreeSamplesFemale.json",
        client = "../data/gen3-client"
    output:
        expand(
        "/project2/yangili1/cfbuenabadn/gtex-stm/code/gtex-download/bams/ThreeSamplesFemale/{IndID}.Aligned.sortedByCoord.out.patched.md.bam", 
        IndID = female_samples
        )
    log:
        'logs/download_SmallTest.log' 
    shell:
        """
        {input.client} download-multiple --profile=AnVIL --manifest={input.manifest} --download-path=/project2/yangili1/cfbuenabadn/gtex-stm/code/gtex-download/bams/ThreeSamplesFemale/ --protocol=s3
        """
        
rule DownloadFromGTEX:
    input:
        manifest = "gtex-download/file_manifest/TrainTest_2tissues.json",
        client = "../data/gen3-client"
    output:
        expand(
        "/project2/yangili1/cfbuenabadn/gtex-stm/code/gtex-download/bams/TrainTest_2tissues/{IndID}.Aligned.sortedByCoord.out.patched.md.bam", 
        IndID = gtex_samples
        )
    log:
        'logs/download_TrainTest_2tissues.log' 
    shell:
        """
        {input.client} download-multiple --profile=AnVIL --manifest={input.manifest} --download-path=/project2/yangili1/cfbuenabadn/gtex-stm/code/gtex-download/bams/TrainTest_2tissues/ --protocol=s3
        """

rule GetIndex:
    input:
        bam = "/project2/yangili1/cfbuenabadn/gtex-stm/code/gtex-download/bams/{SampleSet}/{IndID}.Aligned.sortedByCoord.out.patched.md.bam"
    output:
        bai = "/project2/yangili1/cfbuenabadn/gtex-stm/code/gtex-download/bams/{SampleSet}/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.bai"
    log:
        "logs/bam_idx.{SampleSet}.{IndID}.log"
    shell:
        """
        samtools index {input.bam} {output.bai} &> {log}
        """

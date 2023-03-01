rule MakeFileManifestJson:
    input:
        "../data/file-manifest.json",
        '../data/sample.tsv',
        '../data/participant.tsv'
    output:
        "gtex-download/{Tissue}/files/tissue-manifest.json",
    params:
        GetTissueListStr
    log:
        'logs/makefilemanifest/{Tissue}.json'
    shell:
        """
        python scripts/make_tissue_manifest.py --tissue_samples {params} --output {output} &> {log}
        """
        
rule DownloadFromGTEx_Brain_Cortex:
    input:
        manifest = "gtex-download/{Tissue}/files/tissue-manifest.json",
        client = "../data/gen3-client"
    output:
        expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam", 
        IndID = brain_cortex_samples
        )
    log:
        'logs/download_{Tissue}.log' 
    wildcard_constraints:
        Tissue = 'Brain_Cortex'
    shell:
        """
        {input.client} download-multiple --profile=AnVIL --manifest={input.manifest} --download-path=/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{wildcards.Tissue}/bams/ --protocol=s3
        """
        
use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Muscle_Skeletal with:
    output:
        expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam", 
        IndID = muscle_skeletal_samples
        )
    wildcard_constraints:
        Tissue = 'Muscle_Skeletal'
        
use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Whole_Blood with:
    output:
        expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam", 
        IndID = whole_blood_samples
        )
    wildcard_constraints:
        Tissue = 'Whole_Blood'
        
        
use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Liver with:
    output:
        expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam", 
        IndID = liver_samples
        )
    wildcard_constraints:
        Tissue = 'Liver'

#def GetFileManifest(wildcards):
#    manifest = "gtex_download/file_manifest/{SampleSet}.json".format(SampleSet = wildcards.SampleSet)
#    return manifest

#rule DownloadFromGTEX_ThreeSamplesFemale:
#    input:
#        manifest = "gtex-download/file_manifest/ThreeSamplesFemale.json",
#        client = "../data/gen3-client"
#    output:
#        expand(
#        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/bams/{/{IndID}.Aligned.sortedByCoord.out.patched.md.bam", 
#        IndID = female_samples
#        )
#    log:
#        'logs/download_SmallTest.log' 
#    shell:
#        """
#        {input.client} download-multiple --profile=AnVIL --manifest={input.manifest} --download-path=/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/bams/ThreeSamplesFemale/ --protocol=s3
#        """
        
#rule DownloadFromGTEX:
#    input:
#        manifest = "gtex-download/file_manifest/TrainTest_2tissues.json",
#        client = "../data/gen3-client"
#    output:
#        expand(
#        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/bams/TrainTest_2tissues/{IndID}.Aligned.sortedByCoord.out.patched.md.bam", 
#        IndID = gtex_samples
#        )
#    log:
#        'logs/download_TrainTest_2tissues.log' 
#    shell:
#        """
#        {input.client} download-multiple --profile=AnVIL --manifest={input.manifest} --download-path=/project2/mstephens/cfbuenabadn/gtex-#stm/code/gtex-download/bams/TrainTest_2tissues/ --protocol=s3
#        """


rule BamIndex:
    input:
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{Tissue}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam"
    output:
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{Tissue}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.bai"
    log:
        "logs/indexbam/{Tissue}.{IndID}.log"
    shell:
        """
        samtools index {input} -o {output} > {log}
        """

rule GetIndex:
    input:
        bam = "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/bams/{SampleSet}/{IndID}.Aligned.sortedByCoord.out.patched.md.bam"
    output:
        bai = "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/bams/{SampleSet}/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.bai"
    log:
        "logs/bam_idx.{SampleSet}.{IndID}.log"
    shell:
        """
        samtools index {input.bam} {output.bai} &> {log}
        """

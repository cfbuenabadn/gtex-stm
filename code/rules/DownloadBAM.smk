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
        temp(expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam", 
        IndID = brain_cortex_samples
        ))
    log:
        'logs/download_{Tissue}.log' 
    wildcard_constraints:
        Tissue = 'Brain_Cortex'
    resources:
        mem_mb = 24000
    shell:
        """
        ({input.client} download-multiple --profile=AnVIL --manifest={input.manifest} --download-path=/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{wildcards.Tissue}/bams/ --protocol=s3) &> {log}
        """
        
use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Muscle_Skeletal with:
    output:
        temp(expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam", 
        IndID = muscle_skeletal_samples
        ))
    wildcard_constraints:
        Tissue = 'Muscle_Skeletal'
        
use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Whole_Blood with:
    output:
        temp(expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam", 
        IndID = whole_blood_samples
        ))
    wildcard_constraints:
        Tissue = 'Whole_Blood'
        
        
use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Liver with:
    output:
        temp(expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam", 
        IndID = liver_samples
        ))
    wildcard_constraints:
        Tissue = 'Liver'

use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Brain_Hippocampus with:
    output:
        temp(expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        IndID = brain_hippocampus_samples
        ))
    wildcard_constraints:
        Tissue = 'Brain_Hippocampus'

use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Brain_Hypothalamus with:
    output:
        expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        IndID = brain_hypothalamus_samples
        )
    wildcard_constraints:
        Tissue = 'Brain_Hypothalamus'

use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Brain_Cerebellum with:
    output:
        expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        IndID = brain_cerebellum_samples
        )
    wildcard_constraints:
        Tissue = 'Brain_Cerebellum'
        
use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Kidney_Cortex with:
    output:
        expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        IndID = kidney_cortex_samples
        )
    wildcard_constraints:
        Tissue = 'Kidney_Cortex'
        
use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Lung with:
    output:
        temp(expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        IndID = lung_samples
        ))
    wildcard_constraints:
        Tissue = 'Lung'

use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Spleen with:
    output:
        expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        IndID = spleen_samples
        )
    wildcard_constraints:
        Tissue = 'Spleen'
        
use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Heart_Atrial_Appendage with:
    output:
        expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        IndID = heart_atrial_appendage_samples
        )
    wildcard_constraints:
        Tissue = 'Heart_Atrial_Appendage'
        
use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Skin_Not_Sun_Exposed_Suprapubic with:
    output:
        expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        IndID = skin_not_sun_exposed_suprapubic_samples
        )
    wildcard_constraints:
        Tissue = 'Skin_Not_Sun_Exposed_Suprapubic'

rule BamIndex:
    input:
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{Tissue}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam"
    output:
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{Tissue}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.bai"
    log:
        "logs/indexbam/{Tissue}.{IndID}.log"
    shell:
        """
        samtools index {input} > {log}
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

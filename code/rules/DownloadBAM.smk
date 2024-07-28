rule MakeFileManifestJson:
    input:
        "manifests/file-manifest.json",
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

use rule MakeFileManifestJson as MakeFileManifestJson_subsample with:
    output:
        'gtex-download/subset/files/tissue-manifest.json'
    params:
        '.'.join(subsamples)
    log:
        'logs/makefilemanifest/subsamples_.json'
    
        
rule DownloadFromGTEx_Brain_Cortex:
    input:
        manifest = "gtex-download/{Tissue}/files/tissue-manifest.json",
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
        mem_mb = 42000
    shell:
        """
        (./gen3-client download-multiple --profile=AnVIL --manifest={input.manifest} --download-path=/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{wildcards.Tissue}/bams/ --protocol=s3) &> {log}
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
        temp(expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        IndID = brain_hypothalamus_samples
        ))
    wildcard_constraints:
        Tissue = 'Brain_Hypothalamus'

use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Brain_Cerebellum with:
    output:
        temp(expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        IndID = brain_cerebellum_samples
        ))
    wildcard_constraints:
        Tissue = 'Brain_Cerebellum'
        
use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Kidney_Cortex with:
    output:
        temp(expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        IndID = kidney_cortex_samples
        ))
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
        temp(expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        IndID = spleen_samples
        ))
    wildcard_constraints:
        Tissue = 'Spleen'
        
use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Heart_Atrial_Appendage with:
    output:
        temp(expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        IndID = heart_atrial_appendage_samples
        ))
    wildcard_constraints:
        Tissue = 'Heart_Atrial_Appendage'
        
use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Skin_Not_Sun_Exposed_Suprapubic with:
    output:
        temp(expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        IndID = skin_not_sun_exposed_suprapubic_samples
        ))
    wildcard_constraints:
        Tissue = 'Skin_Not_Sun_Exposed_Suprapubic'
        

use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Cells_Cultured_fibroblasts with:
    output:
        expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        IndID = fibroblast_samples
        )
    wildcard_constraints:
        Tissue = 'Cells_Cultured_fibroblasts'
        
use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_LCL with:
    output:
        expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        IndID = LCL_samples
        )
    wildcard_constraints:
        Tissue = 'Cells_EBV-transformed_lymphocytes'
        
use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Brain_Frontal_Cortex_BA9 with:
    output:
        temp(expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        IndID = BA9_samples
        ))
    wildcard_constraints:
        Tissue = 'Brain_Frontal_Cortex_BA9'
        
use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Brain_Anterior_cingulate_cortex_BA24 with:
    output:
        temp(expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        IndID = BA24_samples
        ))
    wildcard_constraints:
        Tissue = 'Brain_Anterior_cingulate_cortex_BA24'
        
use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Brain_Putamen_basal_ganglia with:
    output:
        temp(expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        IndID = putamen_samples
        ))
    wildcard_constraints:
        Tissue = 'Brain_Putamen_basal_ganglia'

use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Brain_Caudate_basal_ganglia with:
    output:
        temp(expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        IndID = caudate_samples
        ))
    wildcard_constraints:
        Tissue = 'Brain_Caudate_basal_ganglia'
        

use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Brain_Cerebellar_Hemisphere with:
    output:
        temp(expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        IndID = cerebellarh_samples
        ))
    wildcard_constraints:
        Tissue = 'Brain_Cerebellar_Hemisphere'



rule BamIndex:
    input:
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{Tissue}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam"
    output:
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{Tissue}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.bai"
    log:
        "logs/indexbam/{Tissue}.{IndID}.log"
    resources:
        mem_mb = 24000
    shell:
        """
        samtools index {input} {output} > {log}
        """

#rule GetIndex:
#    input:
#        bam = "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/bams/{SampleSet}/{IndID}.Aligned.sortedByCoord.out.patched.md.bam"
#    output:
#        bai = "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-#download/bams/{SampleSet}/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.bai"
#    log:
#        "logs/bam_idx.{SampleSet}.{IndID}.log"
#    shell:
#        """
#        samtools index {input.bam} {output.bai} &> {log}
#        """

rule DownloadFromGTEx_TestSamples:
    input:
        manifest = "manifests/test-manifest.json",
    output:
        expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/TestSamples/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        IndID = test_samples
        )
    log:
        'logs/download_test.log' 
    resources:
        mem_mb = 42000
    shell:
        """
        (./gen3-client download-multiple --profile=AnVIL --manifest={input.manifest} --download-path=/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/TestSamples/bams/ --protocol=s3) &> {log}
        """
        
rule DownloadFromGTEx_TestSamples_juncs:
    input:
        manifest = "manifests/test-manifest-juncs.json",
    output:
        expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/TestSamples/juncs/{IndID}.leafcutter.junc.gz",
        IndID = test_samples
        )
    log:
        'logs/download_test.log' 
    resources:
        mem_mb = 42000
    shell:
        """
        (./gen3-client download-multiple --profile=AnVIL --manifest={input.manifest} --download-path=/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/TestSamples/juncs/ --protocol=s3) &> {log}
        """
        

use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Muscle_Skeletal_rMATS with:
    output:
        temp(expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download-for-DS/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam", 
        IndID = muscle_skeletal_samples
        ))
    wildcard_constraints:
        Tissue = 'Muscle_Skeletal'

use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Brain_Frontal_Cortex_BA9_rMATS with:
    output:
        temp(expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download-for-DS/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        IndID = BA9_samples
        ))
    wildcard_constraints:
        Tissue = 'Brain_Frontal_Cortex_BA9'

use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Whole_Blood_rMATS with:
    output:
        temp(expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download-for-DS/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        IndID = whole_blood_samples
        ))
    wildcard_constraints:
        Tissue = 'Whole_Blood'

use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Liver_rMATS with:
    output:
        temp(expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download-for-DS/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        IndID = liver_samples
        ))
    wildcard_constraints:
        Tissue = 'Liver'

rule BamIndex_DS:
    input:
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download-for-DS/{Tissue}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam"
    output:
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download-for-DS/{Tissue}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.bai"
    log:
        "logs/indexbam-for-DS/{Tissue}.{IndID}.log"
    resources:
        mem_mb = 24000
    shell:
        """
        samtools index {input} {output} > {log}
        """














use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Adipose_Subcutaneous with:
    output:
        temp(expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        IndID = adipose_subcutaneous_samples
        ))
    wildcard_constraints:
        Tissue = 'Adipose_Subcutaneous'


use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Breast_Mammary_Tissue with:
    output:
        temp(expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        IndID = Breast_Mammary_Tissue_samples
        ))
    wildcard_constraints:
        Tissue = 'Breast_Mammary_Tissue'

use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Heart_Left_Ventricle with:
    output:
        temp(expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        IndID = Heart_Left_Ventricle_samples
        ))
    wildcard_constraints:
        Tissue = 'Heart_Left_Ventricle'

use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Skin_Sun_Exposed_Lower_leg with:
    output:
        temp(expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        IndID = Skin_Sun_Exposed_Lower_leg_samples
        ))
    wildcard_constraints:
        Tissue = 'Skin_Sun_Exposed_Lower_leg'

use rule DownloadFromGTEx_Brain_Cortex as DownloadFromGTEx_Testis with:
    output:
        temp(expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        IndID = Testis_samples
        ))
    wildcard_constraints:
        Tissue = 'Testis'



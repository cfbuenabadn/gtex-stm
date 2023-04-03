#### Note: run these after loading Midway's R4.1.0 by running module load R/4.1.0

#def GetSetBAM(wildcards):
#    if wildcards.Test == 'female_test':
#        IndID_list = get_all_samples(wildcards.Tissue, True)[0]
#    elif wildcards.Test == 'train':
#        IndID_list = get_all_samples(wildcards.Tissue, True)[1]
#    elif wildcards.Test == 'test':
#        IndID_list = get_all_samples(wildcards.Tissue, True)[2]
#        
#    x = expand(
#        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam", 
#        IndID = IndID_list
#        )
#    
#    return x
    
#def GetSetJunc(wildcards):
#    if wildcards.Test == 'ThreeSamplesFemale':
#        IndID_list = get_all_samples(wildcards.Tissue, True)[0]
#    elif wildcards.Test == 'Train':
#        IndID_list = get_all_samples(wildcards.Tissue, True)[1]
#    elif wildcards.Test == 'Test':
#        IndID_list = get_all_samples(wildcards.Tissue, True)[2]
#        
#    x = expand(
#        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc", 
#        IndID = IndID_list
#        )
#    
#    return x
    

    
#def GetCondition1(wildcards):
#    if wildcards.Test == 'TrainTest_2tissues':
#        x = '.'.join(brain_train)
#    elif wildcards.Test == 'ThreeSamplesFemale':
#        x = '.'.join(brain_female_three_samples_only)
#    elif wildcards.Test == 'ThreeSamplesMale':
#        x = '.'.join(brain_male_three_samples_only)
#    return x
    
#def GetCondition2(wildcards):
#    if wildcards.Test == 'TrainTest_2tissues':
#        x = '.'.join(muscle_train)
#    elif wildcards.Test == 'ThreeSamplesFemale':
#        x = '.'.join(muscle_female_three_samples_only)
#    elif wildcards.Test == 'ThreeSamplesMale':
#        x = '.'.join(muscle_male_three_samples_only)
#    return x
    

#def GetTissueTestList(wildcards):
#    tissue_samples_tuple = get_all_samples(wildcards.Tissue, True)
#        
#    if wildcards.Test == 'female_test':
#        tissue_samples = tissue_samples_tuple[0]
#    elif wildcards.Test == 'train':
#        tissue_samples = tissue_samples_tuple[1]
#    elif wildcards.Test == 'test':
#        tissue_samples = tissue_samples_tuple[2]
#        
#    tissue_samples = '.'.join(tissue_samples)
    


rule MakeBfiles:
    input:
        GetTissueGroupBAM
    output:
        "gtex-download/{Tissue}/files/rmats_{Group}.b.txt",
    params:
        tissue_list = GetTissueGroupListStr
    log:
        "logs/differential_splicing/{Group}.{Tissue}.log"
    wildcard_constraints:
        Group = 'male_test|male_test2|female_test|female_test2|test|train',
        Tissue = '|'.join(tissue_list)
    shell:
        """
        python scripts/make_rmats_input_files.py --group {wildcards.Group} --tissue {wildcards.Tissue} --tissue_list {params.tissue_list} --output {output} &> {log}
        """

rule RunRMATS:
    input:
        b1 = "gtex-download/{Tissue_1}/files/rmats_{Group}.b.txt",
        b2 = "gtex-download/{Tissue_2}/files/rmats_{Group}.b.txt",
        gtf = "Annotations/gencode.v43.basic.annotation.gtf"
    output:
        expand('DifferentialSplicing/rMATS/{{Tissue_1}}_v_{{Tissue_2}}_{{Group}}/output/{output_file}', output_file = rsem_output_list),
        #temp(directory('/scratch/midway2/cnajar/rmats/{Tissue_1}_v_{Tissue_2}_{Group}'))
    log:
        'logs/differential_splicing/{Tissue_1}_v_{Tissue_2}_{Group}.rmats_run.log',
    wildcard_constraints:
        Group = 'male_test|male_test2|female_test|female_test2|test|train',
        Tissue_1 = '|'.join(tissue_list),
        Tissue_2 = '|'.join(tissue_list)
    resources:
        mem_mb = 62000
    conda:
        "../envs/rmats.yml"
    shell:
        """
        python  {config[rMATS]} --b1 {input.b1} --b2 {input.b2} --gtf {input.gtf} -t single --readLength 76 --nthread 8 --od DifferentialSplicing/rMATS/{wildcards.Tissue_1}_v_{wildcards.Tissue_2}_{wildcards.Group}/output/ --tmp /scratch/midway2/cnajar/rmats/{wildcards.Tissue_1}_v_{wildcards.Tissue_2}_{wildcards.Group} &> {log}
        """

rule RunRMATS_negatives:
    input:
        b1 = "gtex-download/{Tissue}/files/rmats_{Sex}_test.b.txt",
        b2 = "gtex-download/{Tissue}/files/rmats_{Sex}_test2.b.txt",
        gtf = "Annotations/gencode.v43.basic.annotation.gtf"
    output:
        expand('DifferentialSplicing/rMATS/negative_tests/{{Tissue}}_{{Sex}}_test/output/{output_file}', output_file = rsem_output_list),
    log:
        'logs/differential_splicing/{Tissue}_{Sex}_test.rmats_run.log'
    wildcard_constraints:
        Tissue = '|'.join(tissue_list),
        Sex = 'female|male'
    resources:
        mem_mb = 12000
    conda:
        '../envs/rmats.yml'
    shell:
        """
        python  {config[rMATS]} --b1 {input.b1} --b2 {input.b2} --gtf {input.gtf} -t single --readLength 76 --nthread 8 --od DifferentialSplicing/rMATS/{wildcards.Tissue}_{wildcards.Sex}_test/output/ --tmp /scratch/midway2/cnajar/rmats/{wildcards.Tissue}_{wildcards.Sex}_test &> {log}
        """



rule MakeJuncFiles_Brain_Cortex:
    input:
        GetTissueBAM
    output:
        expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc", 
        IndID = brain_cortex_samples
        ) 
    log:
        "logs/leafcutter/input/{Tissue}.juncfiles.log"
    resources:
        mem_mb = 24000
    wildcard_constraints:
        Tissue = 'Brain_Cortex'
    shell:
        """
        echo 'Transforming bam to junc' > {log}
        for bamfile in $(ls gtex-download/{wildcards.Tissue}/bams/*.Aligned.sortedByCoord.out.patched.md.bam); do
            (echo Converting $bamfile to $bamfile.junc) &>> {log}
            (samtools index $bamfile) &>> {log}
            (regtools junctions extract -s 0 -a 8 -m 50 -M 500000 $bamfile -o $bamfile.junc) &>> {log}
        done;
        """
        
    
use rule MakeJuncFiles_Brain_Cortex as MakeJuncFiles_Muscle_Skeletal with:
    output:
        expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc", 
        IndID = muscle_skeletal_samples
        ) 
    wildcard_constraints:
        Tissue = 'Muscle_Skeletal'  
        
        
use rule MakeJuncFiles_Brain_Cortex as MakeJuncFiles_Whole_Blood with:
    output:
        expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc", 
        IndID = whole_blood_samples
        ) 
    wildcard_constraints:
        Tissue = 'Whole_Blood'  
        
use rule MakeJuncFiles_Brain_Cortex as MakeJuncFiles_Liver with:
    output:
        expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc", 
        IndID = liver_samples
        ) 
    wildcard_constraints:
        Tissue = 'Liver'  
       
use rule MakeJuncFiles_Brain_Cortex as MakeJuncFiles_Brain_Hippocampus with:
    output:
        expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/{{Tissue}}/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc",
        IndID = brain_hippocampus_samples
        )
    wildcard_constraints:
        Tissue = 'Brain_Hippocampus' 
        
rule MakeLeafcutterInputJuncFiles:
    input:
        Tissue_1_Junc = GetTissueJunc,
        Tissue_2_Junc =GetTissueJunc2
        #Tissue_1_Junc = GetTissueGroupJunc,
        #Tissue_2_Junc =GetTissueGroupJunc2
    output:
        'DifferentialSplicing/leafcutter/juncfiles/{Tissue}_v_{Tissue_2}/test_juncfiles.txt',
    log:
        "logs/leafcutter/input/{Tissue}_v_{Tissue_2}.test_juncfiles.log"
    wildcard_constraints:
        Tissue = '|'.join(tissue_list),
        Tissue_2 = '|'.join(tissue_list)
    shell:
        """
        echo 'Creating test_juncfiles.txt' > {log}
        for bamfile in {input.Tissue_1_Junc}; do
            (echo $bamfile >> {output}) &>> {log}
        done;
        for bamfile in {input.Tissue_2_Junc}; do
            (echo $bamfile >> {output}) &>> {log}
        done;
        """

rule MakeLeafcutterInputGroupsFile:
    input:
        'DifferentialSplicing/leafcutter/juncfiles/{Tissue}_v_{Tissue_2}/test_juncfiles.txt'
    output:
        groups_file = 'DifferentialSplicing/leafcutter/{Tissue}_v_{Tissue_2}_{Group}/input/groups_file.txt'
    log:
        "logs/leafcutter/input/{Tissue}_v_{Tissue_2}_{Group}.groups_file.log"
    params:
        tissue1_list = GetTissueGroupListStr,
        tissue2_list = GetTissueGroupListStr2
    wildcard_constraints:
        Group = 'male_test|male_test2|female_test|female_test2|test|train',
        Tissue = '|'.join(tissue_list),
        Tissue_2 = '|'.join(tissue_list)
    shell:
        """
        python scripts/make_leafcutter_input_files.py --group {wildcards.Group} --tissue1 {wildcards.Tissue} --tissue2 {wildcards.Tissue_2} --tissue1_list {params.tissue1_list} --tissue2_list {params.tissue2_list} --output {output} &> {log}
        """

rule LeafcutterCluster:
    input:
        'DifferentialSplicing/leafcutter/juncfiles/{Tissue}_v_{Tissue_2}/test_juncfiles.txt'
    output:
        'DifferentialSplicing/leafcutter/juncfiles/{Tissue}_v_{Tissue_2}/leafcutter_perind_numers.counts.gz'
    log:
        "logs/leafcutter/run/{Tissue}_v_{Tissue_2}.cluster.log"
    wildcard_constraints:
        Tissue = '|'.join(tissue_list),
        Tissue_2 = '|'.join(tissue_list)
    resources:
        mem_mb = 24000
    shell:
        """
        python scripts/davidaknowles_leafcutter/scripts/leafcutter_cluster_regtools_py3.py -j {input} -m 50 -o DifferentialSplicing/leafcutter/juncfiles/{wildcards.Tissue}_v_{wildcards.Tissue_2}/leafcutter -l 500000 &> {log}
        """
        
rule AnnotateClusters:
    input:
        gtf = 'scripts/davidaknowles_leafcutter/leafcutter/data/gencode.v31.exons.txt.gz',
        clusters = 'DifferentialSplicing/leafcutter/juncfiles/{Tissue}_v_{Tissue_2}/leafcutter_perind_numers.counts.gz',
    output:
        'DifferentialSplicing/leafcutter/juncfiles/{Tissue}_v_{Tissue_2}/cluster_annotation.tab.gz'
    log:
        'logs/leafcutter/annotate_clusters_{Tissue}_v_{Tissue_2}.log'
    shell:
        """
        python scripts/annotate_clusters.py --gtf {input.gtf} --clusters {input.clusters} --output {output} &> {log}
        """
        
        
rule LeafcutterDS:
    input:
        perind_counts = 'DifferentialSplicing/leafcutter/juncfiles/{Tissue}_v_{Tissue_2}/leafcutter_perind_numers.counts.gz',
        groups_file = 'DifferentialSplicing/leafcutter/{Tissue}_v_{Tissue_2}_{Group}/input/groups_file.txt',
        exons_file = 'scripts/davidaknowles_leafcutter/leafcutter/data/gencode31_exons.txt.gz'
    output:
        "DifferentialSplicing/leafcutter/{Tissue}_v_{Tissue_2}_{Group}/output/leafcutter_effect_sizes.txt",
        "DifferentialSplicing/leafcutter/{Tissue}_v_{Tissue_2}_{Group}/output/leafcutter_cluster_significance.txt",
    wildcard_constraints:
        Group = 'male_test|male_test2|female_test|female_test2|test|train',
        Tissue = '|'.join(tissue_list),
        Tissue_2 = '|'.join(tissue_list)
    resources:
        mem_mb = 24000
    log:
        "logs/leafcutter/run/{Tissue}_v_{Tissue_2}_{Group}.ds.log"
    shell:
        """
        Rscript scripts/davidaknowles_leafcutter/scripts/leafcutter_ds.R --num_threads 4 {input.perind_counts} {input.groups_file} -i 3 -o DifferentialSplicing/leafcutter/{wildcards.Tissue}_v_{wildcards.Tissue_2}_{wildcards.Group}/output/leafcutter -e {input.exons_file} &> {log}
        """

        
        
        
        
        











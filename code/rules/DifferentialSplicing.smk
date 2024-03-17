#rule UnzipJuncs:
#    input:
#        expand(
#            "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/TestSamples/juncs/{IndID}.leafcutter.junc.gz",
#            IndID = test_samples
#            )
#    output:
#        expand(
#            "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/TestSamples/juncs/{IndID}.leafcutter.junc",
#            IndID = test_samples
#            )
#    shell:
#        """
#        gunzip /project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/TestSamples/juncs/*
#        """

rule MakeJuncsFromBam:
    input:
        bam = "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/TestSamples/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        bai = "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/TestSamples/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.bai"
    output:
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/TestSamples/juncs/{IndID}.leafcutter.junc",
    log:
         "logs/leafcutter/input/{IndID}.juncfiles.log"
    resources:
        mem_mb = 24000
    wildcard_constraints:
        IndID = '|'.join(test_samples)
    shell:
        """
        echo 'Transforming bam to junc' > {log}
        (echo Converting {input.bam} to {output}) &>> {log}
        (regtools junctions extract -s XS -a 8 -m 50 -M 500000 {input.bam} -o {output}) &>> {log}
        """

rule MakeLeafCutterInput:
    input:
        expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/TestSamples/juncs/{IndID}.leafcutter.junc",
        IndID = test_samples
        )
    output:
        expand('DifferentialSplicing/leafcutter/DS_tests/{tissue_test}_{n}Samples/input/test_juncfiles.txt',
        tissue_test = test_list, n = ['3', '5']),
        expand('DifferentialSplicing/leafcutter/DS_tests/{tissue_test}_{n}Samples/input/groups_file.txt', 
        tissue_test = test_list, n = ['3', '5'])
    log:
        "logs/leafcutter/make_inputs.log"
    shell:
        """
        python scripts/make_leafcutter_input_files.py &> {log}
        """

def LeafcutterDS_params(wildcards):
    if wildcards.n == '3':
        return "-i 2 -g 2"
    elif wildcards.tissue_test == 'BA9_1_v_BA9_2':
        return "-i 3 -g 3"
    else:
        return ""

rule LeafcutterCluster:
    input:
        'DifferentialSplicing/leafcutter/DS_tests/{tissue_test}_{n}Samples/input/test_juncfiles.txt'
    output:
        'DifferentialSplicing/leafcutter/DS_tests/{tissue_test}_{n}Samples/input/leafcutter_perind_numers.counts.gz'
    log:
        "logs/leafcutter/run/cluster.{tissue_test}_{n}Samples.log"
    resources:
        mem_mb = 24000
    wildcard_constraints:
        tissue_test = '|'.join(test_list),
        n = '3|5'
    shell:
        """
        python scripts/davidaknowles_leafcutter/scripts/leafcutter_cluster_regtools_py3.py -j {input} -m 50 -o DifferentialSplicing/leafcutter/DS_tests/{wildcards.tissue_test}_{wildcards.n}Samples/input/leafcutter -l 500000 &> {log}
        """
        
rule AnnotateClusters:
    input:
        gtf = 'scripts/davidaknowles_leafcutter/leafcutter/data/gencode.v31.exons.txt.gz',
        clusters = 'DifferentialSplicing/leafcutter/DS_tests/{tissue_test}_{n}Samples/input/leafcutter_perind_numers.counts.gz',
    output:
        'DifferentialSplicing/leafcutter/DS_tests/{tissue_test}_{n}Samples/input/cluster_annotation.tab.gz'
    log:
        'logs/leafcutter/annotate_clusters.{tissue_test}_{n}Samples.log'
    wildcard_constraints:
        tissue_test = '|'.join(test_list),
        n = '3|5'
    shell:
        """
        python scripts/annotate_clusters.py --gtf {input.gtf} --clusters {input.clusters} --output {output} &> {log}
        """
        
rule MakeExonAnnotationForLeafcutter:
    input:
        "Annotations/gencode.v44.primary_assembly.annotation.gtf"
    output:
        "Annotations/gencode44_exons.txt.gz"
    log:
        "logs/leafcutter/exon_annotation.log"
    shell:
        """
        (echo "chr start end strand gene_name" > Annotations/gencode44_exons.txt) &> {log};
        (awk '$3=="exon" {{print $1, $4, $5, $7, $10}}' {input} - | awk -F'["\.]' '{{print $1, $2}}' - >> Annotations/gencode44_exons.txt) &> {log}
        (gzip Annotations/gencode44_exons.txt) &> {log}
        """
        
rule LeafcutterDS:
    input:
        perind_counts = 'DifferentialSplicing/leafcutter/DS_tests/{tissue_test}_{n}Samples/input/leafcutter_perind_numers.counts.gz',
        groups_file = 'DifferentialSplicing/leafcutter/DS_tests/{tissue_test}_{n}Samples/input/groups_file.txt',
        exons_file = 'Annotations/gencode44_exons.txt.gz' #'scripts/davidaknowles_leafcutter/leafcutter/data/gencode31_exons.txt.gz'
    output:
        "DifferentialSplicing/leafcutter/DS_tests/{tissue_test}_{n}Samples/output/leafcutter_effect_sizes.txt",
        "DifferentialSplicing/leafcutter/DS_tests/{tissue_test}_{n}Samples/output/leafcutter_cluster_significance.txt",
    wildcard_constraints:
        tissue_test = '|'.join(test_list),
        n = '3|5'
    params:
        LeafcutterDS_params
    resources:
        mem_mb = 24000
    log:
        "logs/leafcutter/run/{tissue_test}_{n}.ds.log"
    shell:
        """
        Rscript scripts/davidaknowles_leafcutter/scripts/leafcutter_ds.R --num_threads 4 {input.perind_counts} {input.groups_file} {params} -o DifferentialSplicing/leafcutter/DS_tests/{wildcards.tissue_test}_{wildcards.n}Samples/output/leafcutter -e {input.exons_file} &> {log}
        """

rule collect_leafcutter:
    input:
        expand("DifferentialSplicing/leafcutter/DS_tests/{tissue_test}_{n}Samples/output/leafcutter_effect_sizes.txt",
            tissue_test = test_list, n=['3', '5']),
        expand("DifferentialSplicing/leafcutter/DS_tests/{tissue_test}_{n}Samples/output/leafcutter_cluster_significance.txt",
            tissue_test = test_list, n=['3', '5']),
        expand('DifferentialSplicing/leafcutter/DS_tests/{tissue_test}_{n}Samples/input/cluster_annotation.tab.gz',
            tissue_test = test_list, n=['3', '5'])
        
collect1 = ['Brain_Frontal_Cortex_BA9_v_Muscle_Skeletal',
             'Liver_v_Whole_Blood',
             'Lung_v_Skin_Not_Sun_Exposed_Suprapubic']
             
collect2 = ['Brain_Frontal_Cortex_BA9_v_Brain_Putamen_basal_ganglia']
collect3 = ['BA9_1_v_BA9_2']

rule collect_leafcutter_inputs_3_1:
    input:
        expand("DifferentialSplicing/leafcutter/DS_tests/{tissue_test}_3Samples/input/leafcutter_perind_numers.counts.gz",
        tissue_test = collect1)
        
rule collect_leafcutter_inputs_3_2:
    input:
        expand("DifferentialSplicing/leafcutter/DS_tests/{tissue_test}_3Samples/input/leafcutter_perind_numers.counts.gz",
        tissue_test = collect2)
        
rule collect_leafcutter_inputs_3_3:
    input:
        expand("DifferentialSplicing/leafcutter/DS_tests/{tissue_test}_3Samples/input/leafcutter_perind_numers.counts.gz",
        tissue_test = collect3)
        
rule collect_leafcutter_inputs_5_1:
    input:
        expand("DifferentialSplicing/leafcutter/DS_tests/{tissue_test}_5Samples/input/leafcutter_perind_numers.counts.gz",
        tissue_test = collect1)
        
rule collect_leafcutter_inputs_5_2:
    input:
        expand("DifferentialSplicing/leafcutter/DS_tests/{tissue_test}_5Samples/input/leafcutter_perind_numers.counts.gz",
        tissue_test = collect2)
        
rule collect_leafcutter_inputs_5_3:
    input:
        expand("DifferentialSplicing/leafcutter/DS_tests/{tissue_test}_5Samples/input/leafcutter_perind_numers.counts.gz",
        tissue_test = collect3)
        










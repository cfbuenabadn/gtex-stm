rmats_group_list = ['Brain_Putamen_basal_ganglia', 
                    'Muscle_Skeletal', 
                    'Liver', 
                    'Whole_Blood',
                    'Skin_Not_Sun_Exposed_Suprapubic', 
                    'Lung', 
                    'Brain_Frontal_Cortex_BA9', 
                    'BA9_1', 
                    'BA9_2']

rule MakeBfiles:
    input:
        expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/TestSamples/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        IndID = test_samples
        ),
        expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/TestSamples/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.bai",
        IndID = test_samples
        ),
    output:
        expand("DifferentialSplicing/rMATS/bfiles/{group}.{n}b.txt", 
        group = rmats_group_list, n = ['3', '5'])
    log:
        "logs/differential_splicing/rmats_bfiles.log"
    resources:
        mem_mb = 4000
    shell:
        """
        python scripts/make_rmats_input_files.py &> {log}
        """

rule RunRMATS:
    input:
        b1 = "DifferentialSplicing/rMATS/bfiles/{group1}.{n}b.txt",
        b2 = "DifferentialSplicing/rMATS/bfiles/{group2}.{n}b.txt",
        gtf = "Annotations/gencode.v44.primary_assembly.annotation.gtf"
    output:
        expand('DifferentialSplicing/rMATS/{{group1}}_v_{{group2}}_{{n}}Samples/output/{output_file}', output_file = rsem_output_list),
    log:
        'logs/differential_splicing/{group1}_v_{group2}.{n}b.rmats_run.log',
    wildcard_constraints:
        group1 = '|'.join(rmats_group_list),
        group2 = '|'.join(rmats_group_list),
        n = '3|5'
    resources:
        mem_mb = 48000
    conda:
        "../envs/rmats.yml"
    shell:
        """
        python  {config[rMATS]} --b1 {input.b1} --b2 {input.b2} --gtf {input.gtf} -t single --readLength 76 --nthread 8 --od DifferentialSplicing/rMATS/{wildcards.group1}_v_{wildcards.group2}_{wildcards.n}Samples/output/ --tmp /scratch/midway2/cnajar/rmats/{wildcards.group1}_v_{wildcards.group2}_{wildcards.n}Samples &> {log}
        """

rule collect_rmats:
    input:
        expand("DifferentialSplicing/rMATS/Brain_Frontal_Cortex_BA9_v_Muscle_Skeletal_{n}Samples/output/summary.txt", n = ['3', '5']),
        expand("DifferentialSplicing/rMATS/Brain_Frontal_Cortex_BA9_v_Brain_Putamen_basal_ganglia_{n}Samples/output/summary.txt", n = ['3', '5']),
        expand("DifferentialSplicing/rMATS/Liver_v_Whole_Blood_{n}Samples/output/summary.txt", n = ['3', '5']),
        expand("DifferentialSplicing/rMATS/Lung_v_Skin_Not_Sun_Exposed_Suprapubic_{n}Samples/output/summary.txt", n = ['3', '5']),
        expand("DifferentialSplicing/rMATS/BA9_1_v_BA9_2_{n}Samples/output/summary.txt", n = ['3', '5'])
        

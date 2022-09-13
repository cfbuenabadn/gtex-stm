#### Note: run these after loading Midway's R4.1.0 by running module load R/4.1.0

def GetTestBAM(wildcards):
    if wildcards.Test == 'ThreeSamplesFemale':
        IndID_list = female_samples
    elif wildcards.Test == 'ThreeSamplesMale':
        IndID_list = male_samples
    elif wildcards.Test == 'TrainTest_2tissues':
        IndID_list = gtex_samples
        
    x = expand(
        "/project2/yangili1/cfbuenabadn/gtex-stm/code/gtex-download/bams/{Test}/{{IndID}}.Aligned.sortedByCoord.out.patched.md.bam".format(Test = wildcards.Test), 
        IndID = IndID_list
        )
    
    return x
    
def GetTestJunc(wildcards):
    if wildcards.Test == 'ThreeSamplesFemale':
        IndID_list = female_samples_only
    elif wildcards.Test == 'ThreeSamplesMale':
        IndID_list = male_samples_only
    elif wildcards.Test == 'TrainTest_2tissues':
        IndID_list = train_samples
        
    x = expand(
        "/project2/yangili1/cfbuenabadn/gtex-stm/code/gtex-download/bams/{Test}/{{IndID}}.Aligned.sortedByCoord.out.patched.md.bam.junc".format(Test = wildcards.Test), 
        IndID = IndID_list
        )
    
    return x
    
def GetCondition1(wildcards):
    if wildcards.Test == 'TrainTest_2tissues':
        x = '.'.join(brain_train)
    elif wildcards.Test == 'ThreeSamplesFemale':
        x = '.'.join(brain_female_three_samples_only)
    elif wildcards.Test == 'ThreeSamplesMale':
        x = '.'.join(brain_male_three_samples_only)
    return x
    
def GetCondition2(wildcards):
    if wildcards.Test == 'TrainTest_2tissues':
        x = '.'.join(muscle_train)
    elif wildcards.Test == 'ThreeSamplesFemale':
        x = '.'.join(muscle_female_three_samples_only)
    elif wildcards.Test == 'ThreeSamplesMale':
        x = '.'.join(muscle_male_three_samples_only)
    return x


rule MakeBfiles:
    input:
        GetTestBAM
    output:
        "DifferentialSplicing/rMATS/{Test}/b1.txt",
        "DifferentialSplicing/rMATS/{Test}/b2.txt",
    params:
        condition1_list = GetCondition1,
        condition2_list = GetCondition2
    log:
        "logs/differential_splicing/{Test}.log"
    shell:
        """
        python scripts/make_ds_input_files.py --test {wildcards.Test} --ds_tool rMATS --condition1 Brain_Cortex --condition2 Muscle_Skeletal --condition1_list {params.condition1_list} --condition2_list {params.condition2_list} &> {log}
        """

rMATS_files = ['{ASE}.MATS.JCEC.txt', '{ASE}.MATS.JC.txt', 'fromGTF.{ASE}.txt', 'fromGTF.novelJunction.{ASE}.txt',
                  'fromGTF.novelSpliceSite.{ASE}.txt', 'JCEC.raw.input.{ASE}.txt', 'JC.raw.input.{ASE}.txt']
splicing_types = ['A3SS', 'A5SS', 'MXE', 'RI', 'SE']
    
output_list = ['summary.txt']
for rf in rMATS_files:
    for ASE in splicing_types:
        new_file = rf.format(ASE=ASE)
        output_list.append(new_file)
            

rule RunRMATS:
    input:
        bams = GetTestBAM,
        b1 = "DifferentialSplicing/rMATS/{Test}/b1.txt",
        b2 = "DifferentialSplicing/rMATS/{Test}/b2.txt",
        gtf = "Annotations/gencode.v34.primary_assembly.annotation.gtf"
    output:
        expand('DifferentialSplicing/rMATS/{{Test}}/output/{output_file}', output_file = output_list),
        temp(directory('/scratch/midway2/cnajar/rmats/{Test}'))
    log:
        'logs/differential_splicing/{Test}.rmats_run.log',
    resources:
        mem_mb = 42000
    conda:
        "../envs/rmats.yml"
    shell:
        """
        python  {config[rMATS]} --b1 {input.b1} --b2 {input.b2} --gtf {input.gtf} -t single --readLength 76 --nthread 8 --od DifferentialSplicing/rMATS/{wildcards.Test}/output/ --tmp /scratch/midway2/cnajar/rmats/{wildcards.Test} &> {log}
        """

rule MakeLeafcutterInputJunc_ThreeSamplesFemale:
    input:
        expand(
        "/project2/yangili1/cfbuenabadn/gtex-stm/code/gtex-download/bams/ThreeSamplesFemale/{IndID}.Aligned.sortedByCoord.out.patched.md.bam", 
        IndID = female_samples
        ),
    output:
       expand(
        "/project2/yangili1/cfbuenabadn/gtex-stm/code/gtex-download/bams/ThreeSamplesFemale/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc", 
        IndID = female_samples
        ) 
    log:
        "logs/leafcutter/input/ThreeSamplesFemale.juncfiles.log"
    resources:
        mem_mb = 24000
    shell:
        """
        echo 'Transforming bam to junc' > {log}
        for bamfile in $(ls /project2/yangili1/cfbuenabadn/gtex-stm/code/gtex-download/bams/ThreeSamplesFemale/*.Aligned.sortedByCoord.out.patched.md.bam); do
            (echo Converting $bamfile to $bamfile.junc) &>> {log}
            (samtools index $bamfile) &>> {log}
            (regtools junctions extract -s 0 -a 8 -m 50 -M 500000 $bamfile -o $bamfile.junc) &>> {log}
        done;
        """


rule MakeLeafcutterInputJunc_TrainTest:
    input:
        expand(
        "/project2/yangili1/cfbuenabadn/gtex-stm/code/gtex-download/bams/TrainTest_2tissues/{IndID}.Aligned.sortedByCoord.out.patched.md.bam", 
        IndID = gtex_samples
        ),
    output:
       expand(
        "/project2/yangili1/cfbuenabadn/gtex-stm/code/gtex-download/bams/TrainTest_2tissues/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc", 
        IndID = gtex_samples
        ) 
    log:
        "logs/leafcutter/input/TrainTest_2tissues.juncfiles.log"
    resources:
        mem_mb = 24000
    shell:
        """
        echo 'Transforming bam to junc' > {log}
        for bamfile in $(ls /project2/yangili1/cfbuenabadn/gtex-stm/code/gtex-download/bams/TrainTest_2tissues/*.Aligned.sortedByCoord.out.patched.md.bam); do
            (echo Converting $bamfile to $bamfile.junc) &>> {log}
            (samtools index $bamfile) &>> {log}
            (regtools junctions extract -s 0 -a 8 -m 50 -M 500000 $bamfile -o $bamfile.junc) &>> {log}
        done;
        """
        
rule MakeLeafcutterInputTestJuncFiles:
    input:
        GetTestJunc
    output:
        'DifferentialSplicing/leafcutter/{Test}/input/test_juncfiles.txt',
    log:
        "logs/leafcutter/input/{Test}.test_juncfiles.log"
    shell:
        """
        echo 'Creating test_juncfiles.txt' > {log}
        for bamfile in {input}; do
            (echo $bamfile >> {output}) &>> {log}
        done;
        """

rule MakeLeafcutterInputGroupsFile:
    input:
        bams = GetTestBAM,
    output:
        groups_file = 'DifferentialSplicing/leafcutter/{Test}/input/groups_file.txt'
    log:
        "logs/leafcutter/input/{Test}.groups_file.log"
    params:
        condition1_list = GetCondition1,
        condition2_list = GetCondition2
    shell:
        """
        python scripts/make_ds_input_files.py --test {wildcards.Test} --ds_tool leafcutter --condition1 Brain_Cortex --condition2 Muscle_Skeletal --condition1_list {params.condition1_list} --condition2_list {params.condition2_list} &> {log}
        """

rule LeafcutterCluster:
    input:
        'DifferentialSplicing/leafcutter/{Test}/input/test_juncfiles.txt',
    output:
        "DifferentialSplicing/leafcutter/{Test}/output/leafcutter_perind_numers.counts.gz"
    log:
        "leafcutter/run/{Test}.cluster.log"
    resources:
        mem_mb = 24000
    shell:
        """
        python scripts/davidaknowles_leafcutter/scripts/leafcutter_cluster_regtools_py3.py -j {input} -m 50 -o DifferentialSplicing/leafcutter/{wildcards.Test}/output/leafcutter -l 500000 &> {log}
        """
        
rule LeafcutterDS:
    input:
        perind_counts = "DifferentialSplicing/leafcutter/{Test}/output/leafcutter_perind_numers.counts.gz",
        groups_file = "DifferentialSplicing/leafcutter/{Test}/input/groups_file.txt"
    output:
        "DifferentialSplicing/leafcutter/{Test}/output/leafcutter_effect_sizes.txt",
        "DifferentialSplicing/leafcutter/{Test}/output/leafcutter_cluster_significance.txt",
    resources:
        mem_mb = 24000
    log:
        "logs/leafcutter/run/{Test}.ds.log"
    shell:
        """
        Rscript scripts/davidaknowles_leafcutter/scripts/leafcutter_ds.R --num_threads 4 {input.perind_counts} {input.groups_file} -i 3 -o DifferentialSplicing/leafcutter/{wildcards.Test}/output/leafcutter -e scripts/davidaknowles_leafcutter/leafcutter/data/gencode.v31.exons.txt.gz &> {log}
        """

        
        
        
        
        











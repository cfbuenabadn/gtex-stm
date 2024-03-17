rule GetFastQFromBam:
    input:
        bam = "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/TestSamples/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
        bai = "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/TestSamples/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.bai"
    output:
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/TestSamples/fastq/{IndID}_1.fastq.gz",
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/TestSamples/fastq/{IndID}_2.fastq.gz",
    log:
         "logs/leafcutter/input/{IndID}.fastq.log"
    resources:
        mem_mb = 24000
    wildcard_constraints:
        IndID = '|'.join(test_samples)
    shell:
        """
        python scripts/run_SamToFastq.py {input.bam} -p gtex-download/TestSamples/fastq/{wildcards.IndID} --jar {config[picard]} &> {log}
        """

rule KallistoIndex:
    input:
        "Annotations/gencode.v44.transcripts.fa.gz"
    output:
        "DifferentialSplicing/kallisto/index/gencode.v44.index"
    log:
        "logs/kallisto/index.log"
    threads: 4
    resources:
        mem_mb = 24000
    shell:
        """
        {config[kallisto]} index -i {output} --threads {threads} {input} &> {log}
        """
        
rule KallistoQuant:
    input:
        index = "DifferentialSplicing/kallisto/index/gencode.v44.index",
        fastq1 = "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/TestSamples/fastq/{IndID}_1.fastq.gz",
        fastq2 = "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/TestSamples/fastq/{IndID}_2.fastq.gz"
    output:
        "DifferentialSplicing/kallisto/quant/{IndID}/abundance.tsv",
        "DifferentialSplicing/kallisto/quant/{IndID}/abundance.h5",
    log:
        "logs/kallisto/quant/{IndID}.log"
    resources:
        mem_mb = 12000
    threads: 4
    shell:
        """
        {config[kallisto]} quant -i {input.index} -o DifferentialSplicing/kallisto/quant/{wildcards.IndID}/ -b 100 -t {threads}  {input.fastq1} {input.fastq2} &> {log}
        """
        
rule MakeSleuthInput:
    input:
        expand("DifferentialSplicing/kallisto/quant/{IndID}/abundance.tsv", IndID=test_samples),
        expand("DifferentialSplicing/kallisto/quant/{IndID}/abundance.h5", IndID=test_samples)
    output:
        expand("DifferentialSplicing/kallisto/sleuth/{tissue_test}_{n}Samples/input/groups_file.txt", tissue_test = test_list, n = ['3', '5'])
    log:
        "logs/kallisto/sleuth_input.log"
    shell:
        """
        python scripts/make_sleuth_input_files.py &> {log}
        """
        
rule RunSleuth:
    input:
        "DifferentialSplicing/kallisto/sleuth/{tissue_test}_{n}Samples/input/groups_file.txt"
    output:
        "DifferentialSplicing/kallisto/sleuth/{tissue_test}_{n}Samples/sleuth_results.tab.gz"
    log:
        "logs/kallisto/sleuth_runs/{tissue_test}_{n}.log"
    resources:
        mem_mb = 24000
    shell:
        """
        Rscript scripts/run_sleuth.R {input} {output} %>% {log}
        """

rule CollectKallisto:
    input:
        expand("DifferentialSplicing/kallisto/quant/{IndID}/abundance.h5", IndID=test_samples),
        expand("DifferentialSplicing/kallisto/sleuth/{tissue_test}_{n}Samples/input/groups_file.txt", tissue_test = test_list, n = ['3', '5']),
        expand("DifferentialSplicing/kallisto/sleuth/{tissue_test}_{n}Samples/sleuth_results.tab.gz", tissue_test = test_list, n = ['3', '5']),
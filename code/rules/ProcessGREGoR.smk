rule junc2bed:
    input:
        "gregor_data/sberger/splicing_data/junction/{IndID}.junc"
    output:
        temp("gregor_data/sberger/tmp/{IndID}.bed")
    log:
        "/scratch/midway3/cnajar/logs/gregor/make_junc_bed_{IndID}.log"
    wildcard_constraints:
        IndID = '|'.join(gregor_samples)
    shell:
        """
        (awk '($1 ~ /^chr/) && ($6=="+" || $6=="-") {{ print $1, $2, $3, "{wildcards.IndID}", $5 }}' OFS='\\t' {input} > {output}) &> {log}
        """

chrom_list = ['chr' + str(i) for i in range(1, 23)] + ['chrX']

rule MakeJunctionBed:
    input:
        expand("gregor_data/sberger/tmp/{IndID}.bed", IndID = gregor_samples)
    output:
        "gregor_data/sberger/processed_data/long_tables/{chrom}.junc.bed.gz"
    log:
        "/scratch/midway3/cnajar/logs/gregor/make_junc_bed.{chrom}.log"
    wildcard_constraints:
        chrom = '|'.join(chrom_list)
    shell:
        """
        for junc_file in {input}
        do
          (awk '$1=="{wildcards.chrom}"'  $junc_file >> gregor_data/sberger/processed_data/long_tables/{wildcards.chrom}.junc.bed) &>> {log}
        done
        (gzip gregor_data/sberger/processed_data/long_tables/{wildcards.chrom}.junc.bed) &>> {log}
        """

rule MakeGREGoRJuncCounts:
    input:
        "gregor_data/sberger/processed_data/junc.bed.gz"
    output:
        "gregor_data/sberger/processed_data/junc_counts.bed.gz"
    log:
        "logs/gregor/make_junc_counts.log"
    resources:
        mem_mb = 62000,
    shell:
        """
        python scripts/make_junc_counts_bed.py &> {log}
        """


rule BGZipJuncs:
    input:
        "gregor_data/sberger/splicing_data/junction/{IndID}.junc"
    output:
        bgz = "gregor_data/sberger/processed_data/junction/{IndID}.junc.gz",
        tbi = "gregor_data/sberger/processed_data/junction/{IndID}.junc.gz.tbi"
    log:
        "logs/gregor/bgzip_junc.{IndID}.log"
    resources:
        mem_mb = 12000
    wildcard_constraints:
        IndID = '|'.join(gregor_samples)
    shell:
        """
        (bgzip -c {input} > {output.bgz}) &> {log};
        (tabix -p bed {output.bgz}) &>> {log}
        """

rule TabixGregorBed:
    input:
        "gregor_data/sberger/splicing_data/bedgraph/{IndID}.bed.gz"
    output:
        "gregor_data/sberger/splicing_data/bedgraph/{IndID}.bed.gz.tbi"
    log:
        "logs/gregor/bgzip_bed.{IndID}.log"
    resources:
        mem_mb = 12000
    wildcard_constraints:
        IndID = '|'.join(gregor_samples)
    shell:
        """
        (tabix -p bed {input}) &>> {log}
        """

rule collect_gregor_juncs:
    input:
        expand("gregor_data/sberger/processed_data/long_tables/{chrom}.junc.bed.gz", chrom=chrom_list),
        #expand("gregor_data/sberger/processed_data/junction/{IndID}.junc.gz.tbi", IndID = gregor_samples),
        expand("gregor_data/sberger/splicing_data/bedgraph/{IndID}.bed.gz.tbi", IndID = gregor_samples)

def get_samples_dataset_gregor(wildcards):
    bedgraph_dir = f'/project2/mstephens/cfbuenabadn/gtex-stm/code/gregor_data/{wildcards.dataset}/splicing_data/bedgraph/'
    samples = sorted([x for x in os.listdir(bedgraph_dir) if ((x[-2:]=='gz') or (x[-2:]=='tbi'))])
    samples_files = expand(bedgraph_dir + '{sample}', sample = samples)
    return samples_files

rule MakeGeneCounts_Gregor:
    input:
        get_samples_dataset_gregor,
    output:
        "gregor_data/{dataset}/counts/{gene}.csv.gz"
    resources:
        mem_mb = 24000,
    log:
        "/scratch/midway3/cnajar/logs/counts_gregor/{dataset}.{gene}.log"
    wildcard_constraints:
        gene = '|'.join(list(selected_genes.gene)),
        dataset = 'sberger'
    shell:
        """
        python scripts/prepare_counts_gregor.py {wildcards.gene} {wildcards.dataset} &> {log}
        """

rule collect_gregor:
    input:
        expand('gregor_data/sberger/counts/{gene}.csv.gz', gene = list(selected_genes.gene))






#####################



rule create_snmf_sberger_bash_scripts:
    output:
        expand("ebpmf_models/sberger_models/bash_scripts/{n}.sh", n=range(0,2001))
    resources:
        mem_mb = 8000
    log:
        'logs/create_bash_sberger.log'
    shell:
        """
        python scripts/create_bash_scripts_sberger.py &> {log}
        """





rule snmf_gregor_chunk:
    input:
        bashscript = "ebpmf_models/sberger_models/bash_scripts/{n}.sh",
    wildcard_constraints:
        n = '|'.join([str(x) for x in range(0,2001)]),
    output:
        directory("ebpmf_models/sberger_models/RDS/{n}")
    resources:
        mem_mb = lambda wildcards, attempt: 32000 if int(attempt) == 1 else 58000
    log:
        "logs/snmf_sberger/{n}.log"
    shell:
        """
        bash {input.bashscript} &> {log}
        """

rule collect_gregor_snmf:
    input:
        expand("ebpmf_models/sberger_models/RDS/{n}", n = range(0, 2001))
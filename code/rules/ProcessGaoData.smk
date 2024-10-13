rule TabixGaoFiles:
    input:
        expand('gao_data/bed_files/{{file_name}}.{suffix}.bed.gz', suffix = ['plus', 'minus', 'plus.second_read', 'minus.second_read'])
    output:
        expand('gao_data/bed_files/{{file_name}}.{suffix}.bed.gz.tbi', suffix = ['plus', 'minus', 
               'plus.second_read', 'minus.second_read'])
    log:
        'logs/gao_data/{file_name}.tabix_bed.log'
    resources:
        mem_mb = 12000
    wildcard_constraints:
        file_name = '|'.join(gao_complete_samples)
    shell:
        """
        for I in {input}
        do
            (tabix -p bed ${{I}}) &> {log}
        done
        """

rule collect_tabix_gao:
    input:
        expand('gao_data/bed_files/{file_name}.{suffix}.bed.gz.tbi', file_name = gao_complete_samples,
            suffix = ['plus', 'minus', 
               'plus.second_read', 'minus.second_read'])

rule PrepareCountsGao:
    #input:
    #    expand('gao_data/bed_files/{file_name}', file_name = gao_files),
    #    expand('gao_data/bed_files/{file_name}.tbi', file_name = gao_files)
    output:
        'gao_data/counts/{gene}.csv.gz'
    log:
        'logs/gao_data/counts/{gene}.log'
    resources:
        mem_mb = 12000
    wildcard_constraints:
        gene = '|'.join(gao_genes)
    shell:
        """
        python scripts/prepare_counts_gao.py {wildcards.gene} &> {log}
        """

rule collectGaoCounts:
    input:
        expand('gao_data/counts/{gene_name}.csv.gz', gene_name = gao_genes[15000:])



rule TabixJuncGaoFiles:
    input:
        'gao_data/junc_files/{file_name}.junc'
    output:
        bgz = 'gao_data/junc_files/{file_name}.junc.gz',
        tbi = 'gao_data/junc_files/{file_name}.junc.gz.tbi'
    log:
        'logs/gao_data/junc_file/{file_name}.gzip.log'
    resources:
        mem_mb = 12000
    wildcard_constraints:
        file_name = '|'.join(gao_complete_samples)
    shell:
        """
        (bgzip -c {input} > {output.bgz}) &> {log};
        (tabix -p bed {output.bgz}) &>> {log}
        """

rule collectJuncTabix:
    input:
        expand('gao_data/junc_files/{file_name}.junc.gz.tbi', file_name = gao_complete_samples)




rule create_snmf_gao_bash_scripts:
    output:
        expand("ebpmf_models/gao_models/bash_scripts/{n}.sh", n=range(0,2001))
    resources:
        mem_mb = 8000
    log:
        'logs/create_bash_gao.log'
    shell:
        """
        python scripts/create_bash_scripts_gao.py &> {log}
        """


rule snmf_gao_chunk:
    input:
        bashscript = "ebpmf_models/gao_models/bash_scripts/{n}.sh",
    wildcard_constraints:
        n = '|'.join([str(x) for x in range(0,2001)]),
    output:
        directory("ebpmf_models/gao_models/RDS/{n}")
    resources:
        mem_mb = lambda wildcards, attempt: 32000 if int(attempt) == 1 else 58000
    log:
        "logs/snmf_gao/{n}.log"
    shell:
        """
        bash {input.bashscript} &> {log}
        """

rule collect_gao_snmf:
    input:
        expand("ebpmf_models/gao_models/RDS/{n}", n = range(0, 2001))





rule create_snmf_gao_bash_scripts_batch1:
    output:
        expand("ebpmf_models/gao_models/bash_scripts/batch1.{n}.sh", n=range(0,1001))
    resources:
        mem_mb = 8000
    log:
        'logs/create_bash_gao_batch1.log'
    shell:
        """
        python scripts/create_bash_scripts_gao_batch1.py &> {log}
        """


rule snmf_gao_chunk_batch1:
    input:
        bashscript = "ebpmf_models/gao_models/bash_scripts/batch1.{n}.sh",
    wildcard_constraints:
        n = '|'.join([str(x) for x in range(0,1001)]),
    output:
        directory("ebpmf_models/gao_models/RDS_batch1/{n}")
    resources:
        mem_mb = lambda wildcards, attempt: 32000 if int(attempt) == 1 else 58000
    log:
        "logs/snmf_gao.batch1/{n}.log"
    shell:
        """
        bash {input.bashscript} &> {log}
        """

rule collect_gao_snmf_batch1:
    input:
        expand("ebpmf_models/gao_models/RDS_batch1/{n}", n = range(0, 1001))


rule get_introns_for_ad_plot:
    output:
        'ebpmf_models/gao_models/tables/AD_IR_plots_1e3.bed.gz',
        #'ebpmf_models/gao_models/tables/AD_IR_plots_1e2.bed.gz',
        #'ebpmf_models/gao_models/tables/AD_IR_plots_5e2.bed.gz'
    resources:
        mem_mb = 12000
    log:
        'logs/print_ad_introns_bed.log'
    shell:
        """
        python prepare_ad_introns_plot.py &> {log}
        """


rule get_introns_for_second_ad_plot:
    output:
        'ebpmf_models/gao_models/tables/second_AD_IR_plots_1e3.bed.gz',
        'ebpmf_models/gao_models/tables/second_AD_IR_plots_1e2.bed.gz',
        'ebpmf_models/gao_models/tables/second_AD_IR_plots_5e2.bed.gz'
    resources:
        mem_mb = 12000
    log:
        'logs/print_second_ad_introns_bed.log'
    shell:
        """
        python scripts/ad_intron_plots_second.py &> {log}
        """


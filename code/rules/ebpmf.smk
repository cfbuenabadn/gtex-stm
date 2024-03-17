rule Run_ebpmf_multiple_models:
    input:
        expand("coverage/counts/{Tissue}/{{gene}}.csv.gz", Tissue = tissue_list)
    output:
        'plots/structure/{gene}.K2.train_ini_fasttopics.png',
        'plots/structure/{gene}.K2.train_ini_sgom.png',
        'plots/structure/{gene}.K3.train_ini_fasttopics.png',
        'plots/structure/{gene}.K3.train_ini_sgom.png',
        #'plots/structure/{gene}.K5.train_ini_fasttopics.png',
        #'plots/structure/{gene}.K5.train_ini_sgom.png',
        #'plots/structure/{gene}.K5.train_log_ini_fasttopics.png',
        #'plots/structure/{gene}.K5.train_log_ini_sgom.png',
        'plots/structure/{gene}.K2.test_ini_fasttopics.png',
        'plots/structure/{gene}.K2.test_ini_sgom.png',
        'plots/structure/{gene}.K3.test_ini_fasttopics.png',
        'plots/structure/{gene}.K3.test_ini_sgom.png',
        #'plots/structure/{gene}.K5.test_ini_fasttopics.png',
        #'plots/structure/{gene}.K5.test_ini_sgom.png',
        #'plots/structure/{gene}.K5.test_log_ini_fasttopics.png',
        #'plots/structure/{gene}.K5.test_log_ini_sgom.png',
        'ebpmf_models/RDS/{gene}.ebpmf.rds'
    log:
        'logs/ebpmf/{gene}.log'
    resources:
        mem_mb = 32000,
    shell:
        """
        {config[Rscript]} scripts/run_ebpmf_models.R {wildcards.gene} {output} &> {log}
        """

def GetLogParams(wildcards):
    if wildcards.log == '_log2':
        return "TRUE"
    elif wildcards.log == '':
        return "FALSE"
    else:
        raise Exception('Error in log wildcards')

rule ebpmf_run:
    input:
        "coverage/count_tables/{gene}.tab.gz"
    output:
        'ebpmf_models/RDS/{gene}.K{K}{log}.ebpmf.rds',
        'ebpmf_models/tables/{gene}/ebpmf_K{K}{log}/EF.tab.gz',
        'ebpmf_models/tables/{gene}/ebpmf_K{K}{log}/EF_smooth.tab.gz',
        'ebpmf_models/tables/{gene}/ebpmf_K{K}{log}/EL.tab.gz',
        'ebpmf_models/tables/{gene}/ebpmf_K{K}{log}/EL_test.tab.gz',
        'ebpmf_models/tables/{gene}/ebpmf_sgom_K{K}{log}/EF.tab.gz',
        'ebpmf_models/tables/{gene}/ebpmf_sgom_K{K}{log}/EF_smooth.tab.gz',
        'ebpmf_models/tables/{gene}/ebpmf_sgom_K{K}{log}/EL.tab.gz',
        'ebpmf_models/tables/{gene}/ebpmf_sgom_K{K}{log}/EL_test.tab.gz',
    log:
        'logs/ebpmf_run/{gene}.K{K}{log}.log'
    resources:
        mem_mb = 32000,
    wildcard_constraints:
        gene = '|'.join(list(selected_genes.gene)),
        K = '2|3|4|5|6|7|8|9|10',
        log = '|_log2'
    params:
        GetLogParams
    shell:
        """
        {config[Rscript]} scripts/run_ebpmf.R {wildcards.gene} {wildcards.K} {params} {output} &> {log}
        """
        
rule ebpmf_run_subset:
    input:
        "coverage/count_tables/{gene}.tab.gz"
    output:
        'ebpmf_models/RDS_{subset}/{gene}.K{K}{log}.ebpmf.rds',
        'ebpmf_models/tables_{subset}/{gene}/ebpmf_K{K}{log}/EF.tab.gz',
        'ebpmf_models/tables_{subset}/{gene}/ebpmf_K{K}{log}/EF_smooth.tab.gz',
        'ebpmf_models/tables_{subset}/{gene}/ebpmf_K{K}{log}/EL.tab.gz',
        'ebpmf_models/tables_{subset}/{gene}/ebpmf_K{K}{log}/EL_test.tab.gz',
        'ebpmf_models/tables_{subset}/{gene}/ebpmf_sgom_K{K}{log}/EF.tab.gz',
        'ebpmf_models/tables_{subset}/{gene}/ebpmf_sgom_K{K}{log}/EF_smooth.tab.gz',
        'ebpmf_models/tables_{subset}/{gene}/ebpmf_sgom_K{K}{log}/EL.tab.gz',
        'ebpmf_models/tables_{subset}/{gene}/ebpmf_sgom_K{K}{log}/EL_test.tab.gz',
    log:
        'logs/ebpmf_run/{gene}_{subset}.K{K}{log}.log'
    resources:
        mem_mb = 32000,
    wildcard_constraints:
        gene = '|'.join(list(selected_genes.gene)),
        K = '2|3|4|5|6|7|8|9|10',
        log = '|_log2',
        subset = 'brain|no_brain'
    params:
        GetLogParams
    shell:
        """
        {config[Rscript]} scripts/run_ebpmf_{wildcards.subset}.R {wildcards.gene} {wildcards.K} {params} {output} &> {log}
        """


rule ebpmf_lm_general:
    input:
        'ebpmf_models/RDS/{gene}.K{K}{log}.ebpmf.rds'
    output:
        'ebpmf_models/tables/{gene}/lm_K{K}{log}.tab.gz',
    resources:
        mem_mb = 12000,
    wildcard_constraints:
        gene = '|'.join(list(selected_genes.gene)),
        log = '|_log2',
        K = '2|3|4|5|6|7|8|9|10'
    log:
        'logs/lm/{gene}_K{K}{log}.log'
    shell:
        """
        {config[Rscript]} scripts/lm_general.R {input} {wildcards.K} {output} &> {log}
        """

rule ebpmf_prueba:
    input:
        expand('ebpmf_models/tables/{gene}/lm{log}.tab.gz', gene = prueba_genes, log=['', '_log2']),
        
rule collect_ebpmf:
    input:
        expand('ebpmf_models/tables/{gene}/lm_K{K}.tab.gz', K = ['2', '3', '4', '5'], gene = prueba_genes + list(selected_genes.gene))
    


rule ebpmf_run_filtered:
    input:
        "coverage/counts_filtered/{gene}.csv.gz"
    output:
        'ebpmf_models/filtered/RDS/{gene}.rds',
    log:
        '/scratch/midway3/cnajar/logs/ebpmf_run/{gene}.log'
    resources:
        mem_mb = 62000,
    wildcard_constraints:
        gene = '|'.join(list(selected_genes.gene)),
    shell:
        """
        {config[Rscript]} scripts/run_ebpmf_filtered.R {wildcards.gene} &> {log}
        """

rule collect_ebpmf_filtered:
    input:
        expand('ebpmf_models/filtered/RDS/{gene}.rds', gene = list(selected_genes.gene)[17000:])

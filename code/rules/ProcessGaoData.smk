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

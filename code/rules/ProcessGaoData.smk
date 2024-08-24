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
#    run:
#        from itertools import cycle
#        from contextlib import redirect_stdout
#        with open(log[0], "w") as log_file:
#            log_file('see if error is in running python')
#
#        with redirect_stdout(log_file):
#            snmf_bashscript_pairs = zip(gao_genes, cycle(output))
#            for gene, out_f in snmf_bashscript_pairs:
#                n = out_f.split('/')[-1].removesuffix('.sh')
#                with open(out_f, 'a') as f:
#                    _ = f.write(f'Rscript scripts/run_ebpmf_ad.R {gene} {n}\n')
#                    lines = ['if [ $? -ne 0 ]; then\n', f'echo "sNMF failed for gene {gene}, passing to the next command."\n', 
#                             'else\n', 'echo "sNMF succeeded for gene {gene}."\n', 'fi\n']
#                    _ = f.writelines(lines)
#    
#            # If there are more output files than accession numbers, the extra
#            # output files won't get made in the previous loop and snakemake will
#            # complain of missing output files. as a fail safe, let's append to
#            # each file in output, in effect making an empty file if a file wasn't
#            # made in the for loop above
#            for f in output:
#                open(f, 'a').close()


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
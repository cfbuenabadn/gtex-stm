rule collect_gao_EL:
    output:
        expand('ebpmf_models/gao_models/tables/snmf{k}_EL.tab.gz', k=[3, 5])
    log:
         'logs/collect_gao_EL.log'
    resources:
        mem_mb = 12000
    shell:
        """
        python scripts/collect_ad_EL.py &> {log}
        """

rule collect_gregor_EL:
    output:
        expand('ebpmf_models/sberger_models/tables/{seq}_snmf{k}_EL.tab.gz', k=[3, 5], seq=['seqmatic', 'invitae'])
    log:
         'logs/collect_sberger_EL.log'
    resources:
        mem_mb = 12000
    shell:
        """
        python scripts/collect_sberger_EL.py &> {log}
        """

rule get_snmf_isoforms_gregor_chunks:
    output:
        'ebpmf_models/sberger_models/tables/tmp/seqmatic.snmf_3.{chunk}.gtf',
        'ebpmf_models/sberger_models/tables/tmp/seqmatic.snmf_5.{chunk}.gtf',
        'ebpmf_models/sberger_models/tables/tmp/invitae.snmf_3.{chunk}.gtf',
        'ebpmf_models/sberger_models/tables/tmp/invitae.snmf_5.{chunk}.gtf',
    resources:
        mem_mb = 12000
    log:
        '/scratch/midway3/cnajar/logs/ebpmf_run/collect_gregor_isoforms.{chunk}.log'
    wildcard_constraints:
        chunk = '|'.join([str(x) for x in range(0,2001)]),
    shell:
        """
        python scripts/collect_gregor_isoforms.py {wildcards.chunk} &> {log}
        """

rule collect_gregor_isoforms:
    input:
        expand('ebpmf_models/sberger_models/tables/tmp/{{dataset}}.snmf_{{K}}.{chunk}.gtf', chunk = range(0,2001))
    output:
        'ebpmf_models/sberger_models/tables/{dataset}.snmf_{K}.gtf'
    wildcard_constraints:
        dataset = 'seqmatic|invitae',
        K = '3|5'
    resources:
        mem_mb = 12000
    shell:
        """
        for I in {input}
        do
           cat ${{I}} >> {output}
        done
        """

rule collect_gregor_all_isoforms:
    input:
        expand('ebpmf_models/sberger_models/tables/{dataset}.snmf_{K}.gtf', dataset = ['seqmatic', 'invitae'], K=[3, 5])










rule get_gregor_isoform_exons:
    input:
        'ebpmf_models/sberger_models/tables/{dataset}.snmf_{K}.gtf'
    output:
        tmp = temp('ebpmf_models/sberger_models/tables/{dataset}.snmf_{K}.exons.bed.gz'),
        bed = 'ebpmf_models/sberger_models/tables/{dataset}.snmf_{K}.exons.sorted.bed.gz'
    wildcard_constraints:
        dataset = 'seqmatic|invitae',
        K = '3|5'
    log:
        '/scratch/midway3/cnajar/logs/ebpmf_run/sberger_models_get_isoforms_{dataset}.snmf_{K}.log'
    shell:
        """
        python scripts/gtf2bed.py {input} {output.tmp} &> {log};
        (zcat {output.tmp} | bedtools sort -i - | bgzip -c > {output.bed}) &> {log}
        """

rule get_merged_gregor_isoforms:
    input:
        'ebpmf_models/sberger_models/tables/{dataset}.snmf_{K}.exons.sorted.bed.gz'
    output:
        tmp = temp('ebpmf_models/sberger_models/tables/{dataset}.snmf_{K}.merged_isoforms.exons.bed'),
        bed = 'ebpmf_models/sberger_models/tables/{dataset}.snmf_{K}.merged_isoforms.exons.sorted.bed.gz'
    wildcard_constraints:
        dataset = 'seqmatic|invitae',
        K = '3|5'
    resources:
        mem_mb = 12000
    log:
        '/scratch/midway3/cnajar/logs/ebpmf_run/sberger.{dataset}.snmf_{K}.merge_isoforms.log'
    shell:
        """
        python scripts/merge_isoforms.py {input} {output.tmp} &> {log};
        (bedtools sort -i {output.tmp} | bgzip -c > {output.bed}) &> {log}
        """

##############

n_slice = 100


rule annotate_snmf_gregor_isoforms:
    input:
        gencode = 'Annotations/gencode.v44.primary_assembly.exons.sorted.bed.gz',
        snmf = 'ebpmf_models/sberger_models/tables/{dataset}.snmf_{K}.merged_isoforms.exons.sorted.bed.gz'
    output:
        temp('ebpmf_models/sberger_models/tables/{dataset}.snmf_{K}_{i}.tab.gz')
    resources:
        mem_mb = 12000
    log:
        '/scratch/midway3/cnajar/logs/sberger_models_annotate_{dataset}.snmf_{K}_{i}.log'
    params:
        n = n_slice
    wildcard_constraints:
        i = '|'.join([str(j) for j in range(1, n_slice + 1)]),
        dataset = 'seqmatic|invitae',
        K = '3|5'
    shell:
        """
        python scripts/annotate_gregor_isoforms.py {input} {wildcards.i} {params} {wildcards.dataset} {wildcards.K} &> {log}
        """

rule collect_annotated_gregor_snmf_isoforms:
    input:
        expand('ebpmf_models/sberger_models/tables/{{dataset}}.snmf_{{K}}_{i}.tab.gz', i = [str(j) for j in range(1, n_slice + 1)])
    output:
        'ebpmf_models/sberger_models/tables/annotated.{dataset}.snmf_{K}.tab.gz'
    resources:
        mem_mb = 24000
    wildcard_constraints:
        dataset = 'seqmatic|invitae',
        K = '3|5'
    shell:
        """
        FIRST=true
        for I in {input}
        do
            if [ "$FIRST" = true ]; then
                zcat "${{I}}" >> ebpmf_models/sberger_models/tables/annotated.{wildcards.dataset}.snmf_{wildcards.K}.tab
                FIRST=false
            else
                zcat "${{I}}" | tail -n +2 >> ebpmf_models/sberger_models/tables/annotated.{wildcards.dataset}.snmf_{wildcards.K}.tab
            fi
        done
        gzip ebpmf_models/sberger_models/tables/annotated.{wildcards.dataset}.snmf_{wildcards.K}.tab
        """

rule collect_snmf_gregor_isoform_annotations:
    input:
        expand('ebpmf_models/sberger_models/tables/annotated.{dataset}.snmf_{K}.tab.gz', dataset = ['seqmatic', 'invitae'], K=[3, 5])





############ Gao


rule get_snmf_isoforms_gao_chunks:
    output:
        'ebpmf_models/gao_models/tables/tmp/snmf_3.{chunk}.gtf',
        'ebpmf_models/gao_models/tables/tmp/snmf_5.{chunk}.gtf',
    resources:
        mem_mb = 12000
    log:
        '/scratch/midway3/cnajar/logs/ebpmf_run/collect_gao_isoforms.{chunk}.log'
    wildcard_constraints:
        chunk = '|'.join([str(x) for x in range(0,2001)]),
    shell:
        """
        python scripts/collect_gao_isoforms.py {wildcards.chunk} &> {log}
        """

rule collect_gao_isoforms:
    input:
        expand('ebpmf_models/gao_models/tables/tmp/snmf_{{K}}.{chunk}.gtf', chunk = range(0,2001))
    output:
        'ebpmf_models/gao_models/tables/snmf_{K}.gtf'
    wildcard_constraints:
        K = '3|5'
    resources:
        mem_mb = 12000
    shell:
        """
        for I in {input}
        do
           cat ${{I}} >> {output}
        done
        """

rule collect_gao_all_isoforms:
    input:
        expand('ebpmf_models/gao_models/tables/snmf_{K}.gtf', K=[3, 5])









rule get_gao_isoform_exons:
    input:
        'ebpmf_models/gao_models/tables/snmf_{K}.gtf'
    output:
        tmp = temp('ebpmf_models/gao_models/tables/snmf_{K}.exons.bed.gz'),
        bed = 'ebpmf_models/gao_models/tables/snmf_{K}.exons.sorted.bed.gz'
    wildcard_constraints:
        K = '3|5'
    log:
        '/scratch/midway3/cnajar/logs/ebpmf_run/gao_models_get_isoforms.snmf_{K}.log'
    shell:
        """
        python scripts/gtf2bed.py {input} {output.tmp} &> {log};
        (zcat {output.tmp} | bedtools sort -i - | bgzip -c > {output.bed}) &> {log}
        """

rule get_merged_gao_isoforms:
    input:
        'ebpmf_models/gao_models/tables/snmf_{K}.exons.sorted.bed.gz'
    output:
        tmp = temp('ebpmf_models/gao_models/tables/snmf_{K}.merged_isoforms.exons.bed'),
        bed = 'ebpmf_models/gao_models/tables/snmf_{K}.merged_isoforms.exons.sorted.bed.gz'
    wildcard_constraints:
        K = '3|5'
    resources:
        mem_mb = 12000
    log:
        '/scratch/midway3/cnajar/logs/ebpmf_run/gao.snmf_{K}.merge_isoforms.log'
    shell:
        """
        python scripts/merge_isoforms.py {input} {output.tmp} &> {log};
        (bedtools sort -i {output.tmp} | bgzip -c > {output.bed}) &> {log}
        """

##############

n_slice = 100


rule annotate_snmf_gao_isoforms:
    input:
        gencode = 'Annotations/gencode.v44.primary_assembly.exons.sorted.bed.gz',
        snmf = 'ebpmf_models/gao_models/tables/snmf_{K}.merged_isoforms.exons.sorted.bed.gz'
    output:
        temp('ebpmf_models/gao_models/tables/snmf_{K}_{i}.tab.gz')
    resources:
        mem_mb = 12000
    log:
        '/scratch/midway3/cnajar/logs/gao_models_annotate.snmf_{K}_{i}.log'
    params:
        n = n_slice
    wildcard_constraints:
        i = '|'.join([str(j) for j in range(1, n_slice + 1)]),
        K = '3|5'
    shell:
        """
        python scripts/annotate_gao_isoforms.py {input} {wildcards.i} {params} {wildcards.K} &> {log}
        """

rule collect_annotated_gao_snmf_isoforms:
    input:
        expand('ebpmf_models/gao_models/tables/snmf_{{K}}_{i}.tab.gz', i = [str(j) for j in range(1, n_slice + 1)])
    output:
        'ebpmf_models/gao_models/tables/annotated.snmf_{K}.tab.gz'
    resources:
        mem_mb = 24000
    wildcard_constraints:
        K = '3|5'
    shell:
        """
        FIRST=true
        for I in {input}
        do
            if [ "$FIRST" = true ]; then
                zcat "${{I}}" >> ebpmf_models/gao_models/tables/annotated.snmf_{wildcards.K}.tab
                FIRST=false
            else
                zcat "${{I}}" | tail -n +2 >> ebpmf_models/gao_models/tables/annotated.snmf_{wildcards.K}.tab
            fi
        done
        gzip ebpmf_models/gao_models/tables/annotated.snmf_{wildcards.K}.tab
        """

rule collect_snmf_gao_isoform_annotations:
    input:
        expand('ebpmf_models/gao_models/tables/annotated.snmf_{K}.tab.gz', K=[3, 5])



########################


rule second_annotation_gao:
    input:
        gencode = 'Annotations/gencode.v44.primary_assembly.exons.sorted.bed.gz',
        snmf = 'ebpmf_models/gao_models/tables/snmf_{K}.merged_isoforms.exons.sorted.bed.gz',
        annotation = 'ebpmf_models/gao_models/tables/annotated.snmf_{K}.tab.gz'
    output:
        'ebpmf_models/gao_models/tables/second_annotation.snmf_{K}.merged_isoforms.tab.gz'
    log:
         'logs/gao_second_annot_{K}.log'
    wildcard_constraints:
        K = '3|5'
    resources:
        mem_mb = 12000
    shell:
        """
        python scripts/second_annotation.py {input} {output} &> {log}
        """


rule second_annotation_gregor:
    input:
        gencode = 'Annotations/gencode.v44.primary_assembly.exons.sorted.bed.gz',
        snmf = 'ebpmf_models/sberger_models/tables/{dataset}.snmf_{K}.merged_isoforms.exons.sorted.bed.gz',
        annotation = 'ebpmf_models/sberger_models/tables/annotated.{dataset}.snmf_{K}.tab.gz'
    output:
        'ebpmf_models/sberger_models/tables/second_annotation.{dataset}.snmf_{K}.merged_isoforms.tab.gz'
    log:
         'logs/gregor_second_annot_{dataset}_{K}.log'
    resources:
        mem_mb = 12000
    wildcard_constraints:
        K = '3|5',
        dataset = 'seqmatic|invitae'
    shell:
        """
        python scripts/second_annotation.py {input} {output} &> {log}
        """


rule GetTranscript_EL_gao:
    input:
        exons_bed = 'ebpmf_models/gao_models/tables/snmf_{K}.merged_isoforms.exons.sorted.bed.gz',
        el_bed = 'ebpmf_models/gao_models/tables/snmf{K}_EL.tab.gz'
    output:
        'ebpmf_models/gao_models/tables/transcript.snmf_{K}.EL.bed.gz'
    resources:
        mem_mb = 12000
    log:
        '/scratch/midway3/cnajar/logs/ebpmf_run_gao/transcript_EL_snmf_{K}.log'
    wildcard_constraints:
        K = '3|5'
    shell:
        """
        python scripts/get_transcript_EL_gregor_and_gao.py {input.exons_bed} {input.el_bed} {output} &> {log}
        """


rule GetTranscript_EL_gregor:
    input:
        exons_bed = 'ebpmf_models/sberger_models/tables/{dataset}.snmf_{K}.merged_isoforms.exons.sorted.bed.gz',
        el_bed = 'ebpmf_models/sberger_models/tables/{dataset}_snmf{K}_EL.tab.gz'
    output:
        'ebpmf_models/sberger_models/tables/transcript.{dataset}.snmf_{K}.EL.bed.gz'
    resources:
        mem_mb = 12000
    log:
        '/scratch/midway3/cnajar/logs/ebpmf_run_gregor/{dataset}transcript_EL_snmf_{K}.log'
    wildcard_constraints:
        K = '3|5'
    shell:
        """
        python scripts/get_transcript_EL_gregor_and_gao.py {input.exons_bed} {input.el_bed} {output} &> {log}
        """

rule collect_transcript_EL_gregor_gao:
    input:
        expand('ebpmf_models/sberger_models/tables/transcript.{dataset}.snmf_{K}.EL.bed.gz', K = [3, 5], 
               dataset = ['seqmatic', 'invitae']),
        expand('ebpmf_models/sberger_models/tables/second_annotation.{dataset}.snmf_{K}.merged_isoforms.tab.gz',K = [3, 5], 
               dataset = ['seqmatic', 'invitae']),
        expand('ebpmf_models/gao_models/tables/transcript.snmf_{K}.EL.bed.gz', K = [3, 5]),
        expand('ebpmf_models/gao_models/tables/second_annotation.snmf_{K}.merged_isoforms.tab.gz', K = [3, 5])
        
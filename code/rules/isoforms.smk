n_slice = 200

rule get_snmf_isoforms_slice:
    output:
        'ebpmf_models/filtered/snmf_{K}/tables/tmp/snmf.{correction_tag}_{i}.gtf'
    resources:
        mem_mb = 12000
    log:
        '/scratch/midway3/cnajar/logs/ebpmf_run/collect_isoforms_{K}_{correction_tag}_{i}.log'
    params:
        n = n_slice
    wildcard_constraints:
        i = '|'.join([str(j) for j in range(1, n_slice + 1)]),
        correction_tag = '3prime_corrected|uncorrected',
        K = '2|3|4|5|10'
    shell:
        """
        python scripts/collect_snmf_isoforms_parallel.py {wildcards.i} {params} {wildcards.correction_tag} {wildcards.K} &> {log}
        """

rule collect_snmf_gtf_slices:
    input:
        expand('ebpmf_models/filtered/snmf_{{K}}/tables/tmp/snmf.{{correction_tag}}_{i}.gtf', i = [str(j) for j in range(1, n_slice + 1)])
    output:
        'ebpmf_models/filtered/snmf_{K}/tables/snmf.{correction_tag}.gtf'
    wildcard_constraints:
        correction_tag = '3prime_corrected|uncorrected',
        K = '2|3|4|5|10'
    shell:
        """
        for I in {input}
        do
           cat ${{I}} >> {output}
        done
        """

rule get_snmf_isoform_exons:
    input:
        'ebpmf_models/filtered/snmf_10/tables/snmf.{correction_tag}.gtf'
    output:
        tmp = temp('ebpmf_models/filtered/snmf_10/tables/snmf.{correction_tag}.exons.bed.gz'),
        bed = 'ebpmf_models/filtered/snmf_10/tables/snmf.{correction_tag}.exons.sorted.bed.gz'
    wildcard_constraints:
        correction_tag = '3prime_corrected|uncorrected'
    log:
        '/scratch/midway3/cnajar/logs/ebpmf_run/gtf_to_bed.{correction_tag}.log'
    shell:
        """
        python scripts/gtf2bed.py {input} {output.tmp} &> {log};
        (zcat {output.tmp} | bedtools sort -i - | bgzip -c > {output.bed}) &> {log}
        """

rule get_gencode_isoform_exons:
    input:
        'Annotations/gencode.v44.primary_assembly.annotation.gtf',
    output:
        tmp = temp('Annotations/gencode.v44.primary_assembly.exons.bed.gz'),
        bed = 'Annotations/gencode.v44.primary_assembly.exons.sorted.bed.gz'
    wildcard_constraints:
        correction_tag = '3prime_corrected|uncorrected'
    log:
        '/scratch/midway3/cnajar/logs/ebpmf_run/gencode_gtf_to_bed.log'
    shell:
        """
        python scripts/get_gencode_exons_and_annotation.py &> {log}
        """

rule annotate_snmf_isoforms:
    input:
        gencode = 'Annotations/gencode.v44.primary_assembly.exons.sorted.bed.gz',
        snmf = 'ebpmf_models/filtered/snmf_10/tables/snmf.{correction_tag}.exons.sorted.bed.gz'
    output:
        temp('ebpmf_models/filtered/snmf_10/tables/tmp/annotated.snmf.{correction_tag}_{i}.tab.gz')
    resources:
        mem_mb = 12000
    log:
        '/scratch/midway3/cnajar/logs/ebpmf_run/annotate_isoforms.{correction_tag}_{i}.log'
    params:
        n = n_slice
    wildcard_constraints:
        i = '|'.join([str(j) for j in range(1, n_slice + 1)]),
        correction_tag = '3prime_corrected|uncorrected'
    shell:
        """
        python scripts/annotate_isoforms.py {input} {wildcards.i} {params} {wildcards.correction_tag} &> {log}
        """

rule collect_annotated_snmf_isoforms:
    input:
        expand('ebpmf_models/filtered/snmf_10/tables/tmp/annotated.snmf.{{correction_tag}}_{i}.tab.gz', i = [str(j) for j in range(1, n_slice + 1)])
    output:
        'ebpmf_models/filtered/snmf_10/tables/annotated.snmf.{correction_tag}.tab.gz'
    resources:
        mem_mb = 24000
    wildcard_constraints:
        correction_tag = '3prime_corrected|uncorrected'
    shell:
        """
        FIRST=true
        for I in {input}
        do
            if [ "$FIRST" = true ]; then
                zcat "${{I}}" >> ebpmf_models/filtered/snmf_10/tables/annotated.snmf.{wildcards.correction_tag}.tab
                FIRST=false
            else
                zcat "${{I}}" | tail -n +2 >> ebpmf_models/filtered/snmf_10/tables/annotated.snmf.{wildcards.correction_tag}.tab
            fi
        done
        gzip ebpmf_models/filtered/snmf_10/tables/annotated.snmf.{wildcards.correction_tag}.tab
        """

rule collect_snmf_isoform_annotations:
    input:
        expand('ebpmf_models/filtered/snmf_10/tables/annotated.snmf.{correction_tag}.tab.gz', correction_tag = ['3prime_corrected', 'uncorrected'])

#for I in {input}
#do
#   zcat ${{I}} >> ebpmf_models/filtered/snmf_10/tables/annotated.snmf.tab
#done



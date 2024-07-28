rule CollectJunctionsFromGTEx:
    input:
        'gtex_tables/GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct.gz'
    output:
        'junctions.tab.gz'
    log:
        'logs/make_junctions_tab.log'
    resources:
        mem_mb = 12000
    shell:
        """
    (zcat {input} | tail -n+3 - | awk -F'\\t' 'BEGIN {{OFS="\\t"}} NR==1 {{
        printf "#chrom\\tstart\\tend\\tgene\\t";
        for(i=3; i<=NF; i++) {{
            printf "%s%s", $i, (i==NF ? "\\n" : "\\t")
        }}
    }} NR>1 {{
        gsub("_", "\\t", $1);
        print $0
    }}' - | sed 's/\\t_/\\t/g') | bgzip /dev/stdin -c > {output}
    tabix -p bed junctions.tab.gz
    """


n_slice = 100

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
        'ebpmf_models/filtered/snmf_{K}/tables/snmf.{correction_tag}.gtf'
    output:
        tmp = temp('ebpmf_models/filtered/snmf_{K}/tables/snmf.{correction_tag}.exons.bed.gz'),
        bed = 'ebpmf_models/filtered/snmf_{K}/tables/snmf.{correction_tag}.exons.sorted.bed.gz'
    wildcard_constraints:
        correction_tag = '3prime_corrected|uncorrected',
        K = '3|5|10'
    log:
        '/scratch/midway3/cnajar/logs/ebpmf_run/snmf_{K}.gtf_to_bed.{correction_tag}.log'
    shell:
        """
        python scripts/gtf2bed.py {input} {output.tmp} &> {log};
        (zcat {output.tmp} | bedtools sort -i - | bgzip -c > {output.bed}) &> {log}
        """

rule get_merged_isoforms:
    input:
        'ebpmf_models/filtered/snmf_{K}/tables/snmf.uncorrected.exons.sorted.bed.gz'
    output:
        tmp = temp('ebpmf_models/filtered/snmf_{K}/tables/snmf.merged_isoforms.exons.bed'),
        bed = 'ebpmf_models/filtered/snmf_{K}/tables/snmf.merged_isoforms.exons.sorted.bed.gz'
    wildcard_constraints:
        K = '3|5|10'
    resources:
        mem_mb = 12000
    log:
        '/scratch/midway3/cnajar/logs/ebpmf_run/snmf_{K}.merge_isoforms.log'
    shell:
        """
        python scripts/merge_isoforms.py {input} {output.tmp} &> {log};
        (bedtools sort -i {output.tmp} | bgzip -c > {output.bed}) &> {log}
        """

rule get_gencode_isoform_exons:
    input:
        'Annotations/gencode.v44.primary_assembly.annotation.gtf',
    output:
        tmp = temp('Annotations/gencode.v44.primary_assembly.exons_annotation.bed.gz'),
        bed = 'Annotations/gencode.v44.primary_assembly.exons.sorted.bed.gz'
    wildcard_constraints:
        correction_tag = '3prime_corrected|uncorrected'
    resources:
        mem_mb = 12000
    log:
        '/scratch/midway3/cnajar/logs/ebpmf_run/gencode_gtf_to_bed.log'
    shell:
        """
        python scripts/get_gencode_exons_and_annotation.py &> {log};
        (bedtools sort -i {output.tmp} | gzip - > {output.bed}) &>> {log}
        """

rule annotate_snmf_isoforms:
    input:
        gencode = 'Annotations/gencode.v44.primary_assembly.exons.sorted.bed.gz',
        snmf = 'ebpmf_models/filtered/snmf_{K}/tables/snmf.{correction_tag}.exons.sorted.bed.gz'
    output:
        temp('ebpmf_models/filtered/snmf_{K}/tables/tmp/annotated.snmf.{correction_tag}_{i}.tab.gz')
    resources:
        mem_mb = 12000
    log:
        '/scratch/midway3/cnajar/logs/ebpmf_run/annotate_isoforms.snmf_{K}.{correction_tag}_{i}.log'
    params:
        n = n_slice
    wildcard_constraints:
        i = '|'.join([str(j) for j in range(1, n_slice + 1)]),
        correction_tag = '3prime_corrected|uncorrected|merged_isoforms',
        K = '3|5|10'
    shell:
        """
        python scripts/annotate_isoforms.py {input} {wildcards.i} {params} {wildcards.correction_tag} {wildcards.K} &> {log}
        """

rule collect_annotated_snmf_isoforms:
    input:
        expand('ebpmf_models/filtered/snmf_{{K}}/tables/tmp/annotated.snmf.{{correction_tag}}_{i}.tab.gz', i = [str(j) for j in range(1, n_slice + 1)])
    output:
        'ebpmf_models/filtered/snmf_{K}/tables/annotated.snmf.{correction_tag}.tab.gz'
    resources:
        mem_mb = 24000
    wildcard_constraints:
        correction_tag = '3prime_corrected|uncorrected|merged_isoforms',
        K = '3|5|10'
    shell:
        """
        FIRST=true
        for I in {input}
        do
            if [ "$FIRST" = true ]; then
                zcat "${{I}}" >> ebpmf_models/filtered/snmf_{wildcards.K}/tables/annotated.snmf.{wildcards.correction_tag}.tab
                FIRST=false
            else
                zcat "${{I}}" | tail -n +2 >> ebpmf_models/filtered/snmf_{wildcards.K}/tables/annotated.snmf.{wildcards.correction_tag}.tab
            fi
        done
        gzip ebpmf_models/filtered/snmf_{wildcards.K}/tables/annotated.snmf.{wildcards.correction_tag}.tab
        """

rule collect_snmf_isoform_annotations:
    input:
        expand('ebpmf_models/filtered/snmf_10/tables/annotated.snmf.{correction_tag}.tab.gz', correction_tag = ['3prime_corrected', 'uncorrected'])

#for I in {input}
#do
#   zcat ${{I}} >> ebpmf_models/filtered/snmf_10/tables/annotated.snmf.tab
#done

rule GetTranscriptEL:
    input:
        exons_bed = 'ebpmf_models/filtered/snmf_{K}/tables/snmf.{correction_tag}.exons.sorted.bed.gz',
        el_bed = 'ebpmf_models/filtered/snmf_{K}/tables/EL.bed.gz'
    output:
        'ebpmf_models/filtered/snmf_{K}/tables/transcript.{correction_tag}.EL.bed.gz'
    resources:
        mem_mb = 12000
    log:
        '/scratch/midway3/cnajar/logs/ebpmf_run/transcript_EL.{correction_tag}_snmf_{K}.log'
    wildcard_constraints:
        correction_tag = '3prime_corrected|uncorrected|merged_isoforms',
        K = '3|5|10'
    shell:
        """
        python scripts/get_transcript_EL.py {input.exons_bed} {input.el_bed} {output} &> {log}
        """

rule collect_transcript_EL:
    input:
        expand('ebpmf_models/filtered/snmf_{K}/tables/transcript.{correction_tag}.EL.bed.gz', K = [3, 5, 10], 
               correction_tag = ['3prime_corrected', 'uncorrected'])





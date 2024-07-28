rule DownloadHg38Ref:
    output:
        "Annotations/gencode.v44.{annot}.annotation.gtf"
    wildcard_constraints:
        annot = 'basic|primary_assembly'
    shell:
        """
        wget -O- https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.{wildcards.annot}.annotation.gtf.gz | zcat > {output}
        """
        
rule MakeTranscriptsBed:
    input:
        "Annotations/gencode.v44.primary_assembly.annotation.gtf"
    output:
        genes = "Annotations/genes.bed",
        transcripts = "Annotations/transcripts.bed",
        genes_appris_principal = "Annotations/genes.appris_principal_transcripts.bed"
    shell:
        """
        awk '$3=="gene"' {input} | awk -F'[\\t"]' '{{print $1, $4, $5, $10, $12, $7}}' OFS='\\t' > {output.genes};
        awk '$3=="transcript"' {input} | grep appris_principal - | awk -F'[\\t"]' '{{print $1, $4, $5, $10, $16, $7, $12}}' OFS='\\t' > {output.transcripts};
        bedtools groupby -g 1,4,5,6 -c 2,3 -o min,max -i {output.transcripts} | awk -v OFS='\t' '{print $1, $5, $6, $2, $3, $4}' > {output.genes_appris_principal}
        """

rule GetAnnotatedGenes:
    input:
        "Annotations/gencode.v44.primary_assembly.annotation.gtf"
    output:
        bed = "Annotations/gencode.v44.primary_assembly.genes.bed.gz",
        tbi = "Annotations/gencode.v44.primary_assembly.genes.bed.gz.tbi"
    log:
        "logs/leafcutter/gene_annotation.primary_assembly.log"
    shell:
        """
        (grep "protein_coding" {input} | awk '($3=="gene" && ($1 ~ /^chr/)) {{print $1, $4, $5, $10, $14, $7}}' OFS="\\t" - |  awk -F'[\\t"\.]' '{{print $1, $2, $3, $5, $9, $11}}' OFS='\\t' - | bedtools sort -i - | awk 'BEGIN {{print "#chrom\\tstart\\tend\\tgene_id\\tgene_type\\tstrand"}} {{print}}'  FS="\\t" OFS="\\t" - >> Annotations/gencode.v44.primary_assembly.genes.bed) &> {log};
        (bgzip Annotations/gencode.v44.primary_assembly.genes.bed) &>> {log};
        (tabix -p bed {output.bed}) &>> {log};
        """

rule GetApprisPrincipalExon:
    input:
        "Annotations/gencode.v44.primary_assembly.annotation.gtf"
    output:
        bed = "Annotations/gencode.v44.primary_assembly.appris_principal_exons.bed.gz",
        tbi = "Annotations/gencode.v44.primary_assembly.appris_principal_exons.bed.gz.tbi"
    log:
        "logs/leafcutter/appris_principal_genes.primary_assembly.log"
    shell:
        """
        (awk '$3=="exon" && ($1 ~ /^chr/)'  OFS='\\t' {input} | grep appris_principal - | awk -F'[\\t"]' '{{print $1, $4, $5, $10, $16, $7, $14, $12, $18}}' OFS='\\t' - | awk 'BEGIN {{OFS="\\t"}} {{gsub(/\..*/, "", $4); gsub(/\..*/, "", $8); print}}' - | bedtools sort -i -  | awk 'BEGIN {{print "#chrom\\tstart\\tend\\tgene_id\\tgene_name\\tstrand\\tgene_type\\ttranscript_id\\ttranscript_type"}} {{print}}' FS="\\t" OFS="\\t" - | bgzip -c - > {output.bed}) &> {log}
        (tabix -p bed {output.bed}) &>> {log};
        """

rule GetApprisPrincipalGenes:
    input:
        "Annotations/gencode.v44.primary_assembly.appris_principal_exons.bed.gz"
    output:
        bed = "Annotations/gencode.v44.primary_assembly.appris_principal_genes.bed.gz",
        tbi = "Annotations/gencode.v44.primary_assembly.appris_principal_genes.bed.gz.tbi"
    log:
        "logs/leafcutter/appris_principal_genes.primary_assembly.log"
    shell:
        """
        (zcat {input} | tail -n+2 - | awk '$7=="protein_coding"' OFS='\\t' - | sort -k4 - | bedtools groupby -g 1,4,5,6 -c 2,3 -o min,max -i - | awk '{{print $1, $5, $6, $2, $3, $4}}' OFS='\\t' - | bedtools sort -i - | awk 'BEGIN {{print "#chrom\\tstart\\tend\\tgene_id\\tgene_name\\tstrand"}} {{print}}' FS='\\t' OFS='\\t' - | bgzip -c - > {output.bed}) &> {log};
        (tabix -p bed {output.bed}) &>> {log};
        """

rule GetAllAnnotatedExons:
    input:
        "Annotations/gencode.v44.primary_assembly.annotation.gtf"
    output:
        bed = "Annotations/gencode.v44.primary_assembly.exons.bed.gz",
        tbi = "Annotations/gencode.v44.primary_assembly.exons.bed.gz.tbi"
    log:
        "logs/leafcutter/exon_annotation.primary_assembly.log"
    shell:
        """
        (awk '$3=="exon" && ($1 ~ /^chr/)'  OFS='\\t' {input} | awk -F'[\\t"]' '{{print $1, $4, $5, $10, $16, $7, $14, $12, $18}}' OFS='\\t' - | awk 'BEGIN {{OFS="\\t"}} {{gsub(/\..*/, "", $4); gsub(/\..*/, "", $8); print}}' - | bedtools sort -i - | awk 'BEGIN {{print "#chrom\\tstart\\tend\\tgene_id\\tgene_name\\tstrand\\tgene_type\\ttranscript_id\\ttranscript_type"}} {{print}}'  FS="\\t" OFS="\\t" - | bgzip -c - > {output.bed}) &> {log};
        (tabix -p bed {output.bed}) &>> {log};
        """



rule Get_expressed_AnnotatedExons:
    input:
        "Annotations/gencode.v44.primary_assembly.annotation.gtf"
    output:
        bed = "Annotations/gencode.v44.primary_assembly.protein_coding_exons.bed.gz",
        tbi = "Annotations/gencode.v44.primary_assembly.protein_coding_exons.bed.gz.tbi"
    log:
        "logs/leafcutter/exon_annotation.primary_assembly.log"
    shell:
        """
        (awk '$3=="exon" && ($1 ~ /^chr/)'  OFS='\\t' {input} | grep protein_coding - | awk -F'[\\t"]' '{{print $1, $4, $5, $10, $16, $7, $14, $12, $18}}' OFS='\\t' - | awk 'BEGIN {{OFS="\\t"}} {{gsub(/\..*/, "", $4); gsub(/\..*/, "", $8); print}}' - | bedtools sort -i - | awk 'BEGIN {{print "#chrom\\tstart\\tend\\tgene_id\\tgene_name\\tstrand\\tgene_type\\ttranscript_id\\ttranscript_type"}} {{print}}'  FS="\\t" OFS="\\t" - | bgzip -c - > {output.bed}) &> {log};
        (tabix -p bed {output.bed}) &>> {log};
        """


        
rule DownloadTranscriptome:
    output:
        "Annotations/gencode.v44.transcripts.fa.gz"
    shell:
        """
        wget -O- https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz > {output}
        """


#rule MakeGenesBed:
#    input:
#        "Annotations/gencode.v34.primary_assembly.annotation.gtf",
#    params:
#        extension = 50
#    output:
#        "Annotations/genes.tab.gz"
#    log:
#        "logs/gene_annotation.tab.log"
#    shell:
#        """
#        awk '$3=="gene"' {input} | awk -F'[\\t"]' '{{print $14, $10, $1":"$4-50"-"$5+50}}' OFS='\\t' | gzip - > {output}
#        #python scripts/gtf2bed.py --gtf {input} --extension {params.extension} --output {output} &> {log}
#        """


rule GetAnnotated_protein_coding_genes:
    input:
        "Annotations/gencode.v44.primary_assembly.annotation.gtf"
    output:
        bed = "Annotations/gencode.v44.primary_assembly.protein_coding.genes.bed.gz",
        tbi = "Annotations/gencode.v44.primary_assembly.protein_coding.genes.bed.gz.tbi"
    log:
        "logs/leafcutter/gene_annotation.primary_assembly.protein_coding.log"
    shell:
        """
        (grep "protein_coding" {input} | awk '($3=="gene" && ($1 ~ /^chr/)) {{print $1, $4, $5, $10, $14, $7}}' OFS="\\t" - |  awk -F'[\\t"\.]' '{{print $1, $2, $3, $5, $9, $11}}' OFS='\\t' - | bedtools sort -i - | awk 'BEGIN {{print "#Chr\\tstart\\tend\\tgene_id\\tgene_name\\tstrand"}} {{print}}'  FS="\\t" OFS="\\t" - >> Annotations/gencode.v44.primary_assembly.protein_coding.genes.bed) &> {log};
        (bgzip Annotations/gencode.v44.primary_assembly.protein_coding.genes.bed) &>> {log};
        (tabix -p bed {output.bed}) &>> {log};
        """


rule GetAnnotated_protein_coding_genes_mane_select:
    input:
        "Annotations/gencode.v44.primary_assembly.annotation.gtf"
    output:
        bed = "Annotations/gencode.v44.primary_assembly.mane_select.genes.bed.gz",
        tbi = "Annotations/gencode.v44.primary_assembly.mane_select.genes.bed.gz.tbi"
    log:
        "logs/leafcutter/gene_annotation.primary_assembly.protein_coding.log"
    shell:
        """
        (grep "protein_coding" {input} | grep MANE_Select - | awk '($3=="gene" && ($1 ~ /^chr/)) {{print $1, $4, $5, $10, $14, $7}}' OFS="\\t" - |  awk -F'[\\t"\.]' '{{print $1, $2, $3, $5, $9, $11}}' OFS='\\t' - | bedtools sort -i - | awk 'BEGIN {{print "#Chr\\tstart\\tend\\tgene_id\\tgene_name\\tstrand"}} {{print}}'  FS="\\t" OFS="\\t" - >> Annotations/gencode.v44.primary_assembly.protein_coding.genes.bed) &> {log};
        (bgzip Annotations/gencode.v44.primary_assembly.protein_coding.genes.bed) &>> {log};
        (tabix -p bed {output.bed}) &>> {log};
        """

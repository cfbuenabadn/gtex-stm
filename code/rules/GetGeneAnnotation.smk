rule DownloadHg38Ref:
    output:
        primary_gtf = "Annotations/gencode.v34.primary_assembly.annotation.gtf"
    shell:
        """
        wget -O- ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.primary_assembly.annotation.gtf.gz | zcat > {output.primary_gtf}
        """

rule DownloadHg38Basic:
    output:
        "Annotations/gencode.v43.basic.annotation.gtf"
    shell:
        """
        wget -O- https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.basic.annotation.gtf.gz | zcat > {output}
        """
        
rule MakeGenesBed:
    input:
        "Annotations/gencode.v34.primary_assembly.annotation.gtf",
    params:
        extension = 50
    output:
        "Annotations/genes.tab.gz"
    log:
        "logs/gene_annotation.tab.log"
    shell:
        """
        awk '$3=="gene"' {input} | awk -F'[\\t"]' '{{print $14, $10, $1":"$4-50"-"$5+50}}' OFS='\\t' | gzip - > {output}
        #python scripts/gtf2bed.py --gtf {input} --extension {params.extension} --output {output} &> {log}
        """

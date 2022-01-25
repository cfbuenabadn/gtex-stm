rule DownloadHg38Ref:
    output:
        primary_gtf = "Annotations/gencode.v34.primary_assembly.annotation.gtf"
    shell:
        """
        wget -O- ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.primary_assembly.annotation.gtf.gz | zcat > {output.primary_gtf}
        """
        
rule MakeGenesBed:
    input:
        "Annotations/gencode.v34.primary_assembly.annotation.gtf",
    params:
        extension = 1000
    output:
        "Annotations/genes.bed"
    log:
        "logs/gene_annotation.bed"
    shell:
        """
        python scripts/gtf2bed.py --gtf {input} --extension {params.extension} --output {output} &> {log}
        """
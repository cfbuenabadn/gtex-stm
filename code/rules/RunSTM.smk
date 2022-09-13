rule RunSGOM:
    input:
       "Counts/{gene}.Counts.csv.gz"
    output:
        "stm_models/{gene}/{gene}.sgom_K{knum}.rds"
    log:
        "logs/{gene}.K{knum}.sgom.log"
    wildcard_constraints:
        gene = '|'.join(genes)
    resources:
        mem_mb = 58000
    shell:
        """
        (Rscript scripts/run_sgom.R {input} {wildcards.knum} {output} || true) &> {log}
        """
        
        
        
        
########################## Older, unused ##############################
        
rule RunSGOM_NoTissue:
    input:
       "Counts/{gene}.Counts.no_{tissue}.csv.gz"
    output:
        "stm_models/{gene}/{gene}.no_{tissue}.sgom_K{knum}.rds"
    wildcard_constraints:
        gene = '|'.join(genes)
    log:
        "logs/{gene}.no_{tissue}.K{knum}.sgom.log"
    resources:
        mem_mb = 58000
    shell:
        """
        (Rscript scripts/run_sgom.R {input} {wildcards.knum} {output} || true) &> {log}
        """
        
rule RunSGOM_OneOfTissue:
    input:
       "Counts/{gene}.Counts.one_of_{tissue}.csv.gz"
    output:
        "stm_models/{gene}/{gene}.one_of_{tissue}.sgom_K{knum}.rds"
    wildcard_constraints:
        gene = '|'.join(genes)
    log:
        "logs/{gene}.one_of_{tissue}.K{knum}.sgom.log"
    resources:
        mem_mb = 58000
    shell:
        """
        (Rscript scripts/run_sgom.R {input} {wildcards.knum} {output} || true) &> {log}
        """
        
rule RunSGOM_MinusOneOfTissue:
    input:
       "Counts/{gene}.Counts.minus_one_of_{tissue}.csv.gz"
    output:
        "stm_models/{gene}/{gene}.minus_one_of_{tissue}.sgom_K{knum}.rds"
    wildcard_constraints:
        gene = '|'.join(genes)
    log:
        "logs/{gene}.one_of_{tissue}.K{knum}.sgom.log"
    resources:
        mem_mb = 58000
    shell:
        """
        (Rscript scripts/run_sgom.R {input} {wildcards.knum} {output} || true) &> {log}
        """
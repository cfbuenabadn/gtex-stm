rule RunSGOM:
    input:
       "Counts/{gene}.Counts.csv.gz"
    output:
        "stm_models/{gene}.sgom_K{knum}.rds"
    log:
        "logs/{gene}.K{knum}.sgom.log"
    resources:
        mem_mb = 58000
    shell:
        """
        (Rscript scripts/run_sgom.R {input} {wildcards.knum} || true) &> {log}
        """
        
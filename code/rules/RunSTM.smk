rule RunSGOM:
    input:
       "Counts/{gene}.Counts.csv.gz"
    output:
        "stm_models/{gene}.sgom_K5.rds"
    log:
        "logs/{gene}.sgom.log"
    shell:
        """
        Rscript scripts/run_sgom.R {input} &> {log}
        """
        
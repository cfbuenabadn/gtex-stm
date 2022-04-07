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
        (Rscript scripts/run_sgom.R {input} {wildcards.knum} {output} || true) &> {log}
        """
        
use rule RunSGOM as RunSGOM_50 with:
    input:
        "Counts_50/{gene}.Counts.csv.gz"
    output:
        "stm_models_50/{gene}.sgom_K{knum}.rds"
    log:
        "logs/{gene}.50_samples.K{knum}.sgom.log"
        
        
use rule RunSGOM as RunSGOM_balanced with:
    input:
        "Counts_balanced/{gene}.Counts.csv.gz"
    output:
        "stm_models_balanced/{gene}.sgom_K{knum}.rds"
    log:
        "logs/{gene}.balanced_samples.K{knum}.sgom.log"
        

use rule RunSGOM as RunSGOM_30_1 with:
    input:
        "Counts_30_1/{gene}.Counts.csv.gz"
    output:
        "stm_models_30_1/{gene}.sgom_K{knum}.rds"
    log:
        "logs/{gene}.30_1_samples.K{knum}.sgom.log"
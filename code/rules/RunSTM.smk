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
        

def GetTissueCountsList(wildcards):
    tissue_list = wildcards.Tissues.split('.')
    return expand("Counts/{tissue}/{gene}.Counts.csv.gz",
                  tissue = tissue_list, gene=wildcards.gene)        


def much_more_mem_after_first_attempt(wildcards, attempt):
    if int(attempt) == 1:
        return 24000
    else:
        return 62000

rule train_ebpmf_2tissues:
    input:
        GetTissueCountsList,
        "config/samples.tsv"
    output:
        "ebpmf_model/train_2tissues/{Tissues}-{gene}.K{K}.ebpmf.rds"
    log:
        "logs/{Tissues}-{gene}.K{K}.ebpmf_train.log"
    wildcard_constraints:
        Tissues = "Brain_Cortex.Muscle_Skeletal|Whole_Blood.Liver|Whole_Blood.Muscle_Skeletal",
        gene = '|'.join(genes),
        K = '2|3|5|10'
    resources:
        mem_mb = much_more_mem_after_first_attempt
    shell:
        """
        Rscript scripts/train_ebpmf.R {wildcards.Tissues} {wildcards.gene} {wildcards.K} {output} &> {log}
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

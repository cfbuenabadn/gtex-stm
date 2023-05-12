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
        "ebpmf_model/train_2tissues/{Tissues}/models/{gene}.K{K}.ebpmf.rds",
        "ebpmf_model/train_2tissues/{Tissues}/plots/{gene}.K{K}.ebpmf.StructurePlot.tissue_id.pdf",
        "ebpmf_model/train_2tissues/{Tissues}/plots/{gene}.K{K}.ebpmf.StructurePlot.sex.pdf",
        "ebpmf_model/train_2tissues/{Tissues}/plots/{gene}.K{K}.ebpmf.StructurePlot.tissue_sex_id.pdf"
        #"ebpmf_model/train_2tissues/{Tissues}/plots/{gene}.K{K}.ebpmf.StructurePlot.png",
    log:
        "logs/ebpmf/train/train_2samples/{Tissues}/{gene}.K{K}.ebpmf_train.log"
    wildcard_constraints:
        Tissues = "Brain_Cortex.Muscle_Skeletal|Whole_Blood.Liver|Whole_Blood.Muscle_Skeletal|Brain_Cortex.Brain_Hypothalamus|Brain_Cortex.Brain_Hippocampus|Brain_Hypothalamus.Brain_Cerebellum",
        gene = '|'.join(genes),
        K = '2|3|5|10'
    resources:
        mem_mb = much_more_mem_after_first_attempt
    shell:
        """
        {config[Rscript]} scripts/train_ebpmf.R {wildcards.Tissues} {wildcards.gene} {wildcards.K} {output} &> {log}
        """
        
rule ebpmf_lm:
    input:
        "ebpmf_model/train_2tissues/{tissue1}.{tissue2}/models/{gene}.K2.ebpmf.rds"
    output:
        rds = "ebpmf_model/train_2tissues_scores/{tissue1}.{tissue2}/{gene}.K2.ebpmf.rds",
    log:
        "logs/ebpmf/train/train_2samples/{tissue1}.{tissue2}/{gene}.K2.lm_test.log"
    wildcard_constraints:
        gene = '|'.join(genes),
        tissue1 = 'Brain_Cortex|Brain_Hypothalamus',
        tissue2 = 'Brain_Hippocampus|Muscle_Skeletal|Brain_Cerebellum|Brain_Hypothalamus'
    resources:
        mem_mb = much_more_mem_after_first_attempt
    shell:
        """
        {config[Rscript]} scripts/ebpmf_test.R {input} {wildcards.tissue1} {wildcards.tissue2} {wildcards.gene} {output} &> {log}
        """
        
        
use rule train_ebpmf_2tissues as train_ebpmf_3tissues with:
    output:
        "ebpmf_model/train_3tissues/{Tissues}-{gene}.K{K}.ebpmf.rds",
        "ebpmf_model/train_3tissues/plots/{Tissues}-{gene}.K{K}.ebpmf.StructurePlot.png",
    wildcard_constraints:
        Tissue = "Brain_Cortex.Brain_Hippocampus.Muscle_Skeletal",
        gene = '|'.join(genes),
        K = '2|3|5|10'
        
use rule train_ebpmf_2tissues as train_ebpmf_5tissues with:
    output:
        "ebpmf_model/train_5tissues/models/{gene}-{Tissues}.K{K}.ebpmf.rds",
        "ebpmf_model/train_5tissues/plots/{gene}-{Tissues}.K{K}.ebpmf.StructurePlot.tissue_id.png",
        "ebpmf_model/train_5tissues/plots/{gene}-{Tissues}.K{K}.ebpmf.StructurePlot.sex.png",
        "ebpmf_model/train_5tissues/plots/{gene}-{Tissues}.K{K}.ebpmf.StructurePlot.tissue_sex_id.png",
    wildcard_constraints:
        Tissues = "Brain_Cortex.Brain_Hippocampus.Muscle_Skeletal.Whole_Blood.Liver",
        gene = '|'.join(genes),
        K = '2|3|4|5|10'
        

        
rule train_ebpmf_10tissues:
    input:
        "Counts/Brain_Cerebellum/{gene}.Counts.csv.gz",
        "Counts/Brain_Cortex/{gene}.Counts.csv.gz",
        "Counts/Brain_Hippocampus/{gene}.Counts.csv.gz",
        "Counts/Heart_Atrial_Appendage/{gene}.Counts.csv.gz",
        "Counts/Kidney_Cortex/{gene}.Counts.csv.gz",
        "Counts/Liver/{gene}.Counts.csv.gz",
        "Counts/Lung/{gene}.Counts.csv.gz",
        "Counts/Muscle_Skeletal/{gene}.Counts.csv.gz",
        "Counts/Skin_Not_Sun_Exposed_Suprapubic/{gene}.Counts.csv.gz",
        "Counts/Whole_Blood/{gene}.Counts.csv.gz",
        "config/samples.tsv"
    output:
        "ebpmf_model/train_10tissues/models/{gene}.K{K}.ebpmf.rds",
        "ebpmf_model/train_10tissues/plots/{gene}.K{K}.ebpmf.StructurePlot.tissue_id.pdf",
        "ebpmf_model/train_10tissues/plots/{gene}.K{K}.ebpmf.StructurePlot.sex.pdf",
        "ebpmf_model/train_10tissues/plots/{gene}.K{K}.ebpmf.StructurePlot.tissue_sex_id.pdf",
    wildcard_constraints:
        gene = '|'.join(genes),
        K = '2|3|4|5|10'
    params:
        Tissues = ".".join(["Brain_Cerebellum", "Brain_Cortex", "Brain_Hippocampus", 
                            "Heart_Atrial_Appendage", "Kidney_Cortex", "Liver", "Lung", 
                            "Muscle_Skeletal", "Skin_Not_Sun_Exposed_Suprapubic", "Whole_Blood"])
    log:
        "logs/ebpmf/train/train_10samples/{gene}.K{K}.log"
    resources:
        mem_mb = much_more_mem_after_first_attempt
    shell:
        """
        {config[Rscript]} scripts/train_ebpmf.R {params.Tissues} {wildcards.gene} {wildcards.K} {output} &> {log}
        """
        
rule plot_ebpmf_model:
    input:
        gtf = '/project2/mstephens/cfbuenabadn/gtex-stm/code/Annotations/gencode.v34.primary_assembly.annotation.gtf',
        rds = "ebpmf_model/train_10tissues/models/{geneName}.K{kfactors}.ebpmf.rds",
    output:
        "ebpmf_model/train_10tissues/plots/{geneName}.K{kfactors}.factors.pdf",
        "ebpmf_model/train_10tissues/plots/{geneName}.K{kfactors}.counts.pdf",
        "ebpmf_model/train_10tissues/plots/{geneName}.K{kfactors}.LF.pdf",
    log:
        "logs/ebpmf/train/train_10samples/{geneName}.K{kfactors}.plots.log"
    wildcard_constraints:
        geneName = '|'.join(genes),
        kfactors = '2|3|4|5|10'
    resources:
        mem_mb = much_more_mem_after_first_attempt
    shell:
        """
        python scripts/ebpmf_model_plots.py --GTF {input.gtf} --RDS {input.rds} --geneName {wildcards.geneName} --kfactors {wildcards.kfactors} {output} &> {log}
        """

use rule plot_ebpmf_model as plot_ebpmf_model_2tissues with:
    input:
        gtf = '/project2/mstephens/cfbuenabadn/gtex-stm/code/Annotations/gencode.v34.primary_assembly.annotation.gtf',
        rds = "ebpmf_model/train_2tissues/{Tissues}/models/{geneName}.K{kfactors}.ebpmf.rds",
    output:
        "ebpmf_model/train_2tissues/{Tissues}/plots/{geneName}.K{kfactors}.factors.pdf",
        "ebpmf_model/train_2tissues/{Tissues}/plots/{geneName}.K{kfactors}.counts.pdf",
        "ebpmf_model/train_2tissues/{Tissues}/plots/{geneName}.K{kfactors}.LF.pdf",
    log:
        "logs/ebpmf/train/train_2samples/{Tissues}/{geneName}.K{kfactors}.plots.log"
    wildcard_constraints:
        geneName = '|'.join(genes),
        kfactors = '2|3|4|5|10',
        Tissues = 'Brain_Cortex.Muscle_Skeletal'


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

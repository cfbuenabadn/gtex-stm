rule MakeTrainingSets:
    input:
        "Counts/{gene}.Counts.csv.gz",
        "../data/sample.tsv",
        "../data/participant.tsv"
    output:
        "Tests/{gene}/Counts/Train.csv.gz",
        "Tests/{gene}/Counts/Test.csv.gz",
    log:
        "logs/{gene}.TrainingSets.log"
    shell:
        """
        python scripts/getTrainingSets.py --gene {wildcards.gene} &> {log}
        """

rule CustomizeCountsTable:
    input:
        "Counts/PKM.Counts.csv.gz",
        "../data/sample.tsv",
        "../data/participant.tsv"
    output:
        "Tests/PKM/Counts/SubsetRegion.csv.gz",
        "Tests/PKM/Counts/SubsetRegion_train.csv.gz",
        "Tests/PKM/Counts/SubsetRegion_test.csv.gz",
        "Tests/PKM/Counts/SubsetTissues.csv.gz",
        "Tests/PKM/Counts/SubsetTissues_train.csv.gz",
        "Tests/PKM/Counts/SubsetTissues_test.csv.gz",
        "Tests/PKM/Counts/SubsetTissues_plus_one.csv.gz",
        "Tests/PKM/Counts/SubsetTissues_train_plus_one.csv.gz"
    log:
        "logs/PKM.Tests.log"
    shell:
        """
        python scripts/customCounts.py &> {log}
        """
        
rule RunSGOM_Train:
    input:
       "Tests/{gene}/Counts/Train.csv.gz"
    output:
        "Tests/{gene}/stm_models/Training.sgom_K{knum}.rds"
    log:
        "logs/{gene}.K{knum}.sgom_training.log"
    resources:
        mem_mb = 58000
    shell:
        """
        (Rscript scripts/run_sgom.R {input} {wildcards.knum} {output} || true) &> {log}
        """

rule RunSGOM_PKM:
    input:
        SubsetRegion = "Tests/PKM/Counts/SubsetRegion.csv.gz",
        SubsetRegion_train = "Tests/PKM/Counts/SubsetRegion_train.csv.gz",
        SubsetTissues = "Tests/PKM/Counts/SubsetTissues.csv.gz",
        SubsetTissues_train = "Tests/PKM/Counts/SubsetTissues_train.csv.gz",
        SubsetTissues_plus_one = "Tests/PKM/Counts/SubsetTissues_plus_one.csv.gz",
        SubsetTissues_train_plus_one = "Tests/PKM/Counts/SubsetTissues_train_plus_one.csv.gz"
    output:
        SubsetRegion = "Tests/PKM/stm_models/SubsetRegion.sgom_K{knum}.rds",
        SubsetRegion_train = "Tests/PKM/stm_models/SubsetRegion_train.sgom_K{knum}.rds",
        SubsetTissues = "Tests/PKM/stm_models/SubsetTissues.sgom_K{knum}.rds",
        SubsetTissues_train = "Tests/PKM/stm_models/SubsetTissues_train.sgom_K{knum}.rds",
        SubsetTissues_plus_one = "Tests/PKM/stm_models/SubsetTissues_plus_one.sgom_K{knum}.rds",
        SubsetTissues_train_plus_one = "Tests/PKM/stm_models/SubsetTissues_train_plus_one.sgom_K{knum}.rds",
    log:
        "logs/PKM_subset_sgom.K{knum}.log"
    resources:
        mem_mb = 58000
    shell:
        """
        (Rscript scripts/run_sgom.R {input.SubsetRegion} {wildcards.knum} {output.SubsetRegion} || true) &>> {log};
        (Rscript scripts/run_sgom.R {input.SubsetRegion_train} {wildcards.knum} {output.SubsetRegion_train} || true) &>> {log};
        (Rscript scripts/run_sgom.R {input.SubsetTissues} {wildcards.knum} {output.SubsetTissues} || true) &>> {log};
        (Rscript scripts/run_sgom.R {input.SubsetTissues_train} {wildcards.knum} {output.SubsetTissues_train} || true) &>> {log};
        (Rscript scripts/run_sgom.R {input.SubsetTissues_plus_one} {wildcards.knum} {output.SubsetTissues_plus_one} || true) &>> {log};
        (Rscript scripts/run_sgom.R {input.SubsetTissues_train_plus_one} {wildcards.knum} {output.SubsetTissues_train_plus_one} || true) &>> {log}
        """



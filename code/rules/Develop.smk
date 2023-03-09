rule gather_counts:
    input:
        expand("Counts/{Tissue}/{gene}.Counts.csv.gz", 
               Tissue = ["Brain_Cortex", "Muscle_Skeletal", "Liver", "Whole_Blood"],
               gene=genes_test)

rule ebpmf_test_run:
    input:
        expand("ebpmf_model/train_2tissues/{Tissues}-{gene}.K{K}.ebpmf.rds",
               Tissues = ["Brain_Cortex.Muscle_Skeletal", "Whole_Blood.Liver", "Whole_Blood.Muscle_Skeletal"],
               gene = genes_test, K = [2, 3, 5, 10])

rule ebpmf_train:
    input:
        expand("ebpmf_model/train_2tissues/{Tissues}-{gene}.K{K}.ebpmf.rds",
               Tissues = ["Brain_Cortex.Muscle_Skeletal"],
               gene = genes, K = [2, 3, 5, 10])


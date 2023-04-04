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
    priority: 1

rule ebpmf_train_light:
    input:
        expand("ebpmf_model/train_3tissues/{Tissues}-{gene}.K{K}.ebpmf.rds",
               Tissues = ["Brain_Cortex.Brain_Hippocampus.Muscle_Skeletal"],
               gene = quick_test_genes, K = [2, 3, 5, 10])

rule ebpmf_train:
    input:
        expand("ebpmf_model/train_2tissues/{Tissues}-{gene}.K{K}.ebpmf.rds",
               Tissues = ["Brain_Cortex.Muscle_Skeletal"],
               gene = genes, K = [2, 3, 5, 10]),
        #expand("ebpmf_model/train_2tissues/{Tissues}-{gene}.K{K}.ebpmf.rds",
        #       Tissues = ["Brain_Cortex.Brain_Hippocampus", "Brain_Hypothalamus.Brain_Cerebellum"],
        #       gene = ['NRXN1', 'NRXN2', 'NRXN3'], K = [2, 3, 5, 10]),

rule KL_test:
    input:
        expand("ebpmf_model/train_2tissues_scores/Brain_Cortex.Muscle_Skeletal-{gene}.K2.ebpmf.rds", gene=genes),
        expand("ebpmf_model/train_2tissues_scores/Brain_Cortex.Brain_Hippocampus-{gene}.K2.ebpmf.rds", gene=['NRXN1', 'NRXN2', 'NRXN3']),
        expand("ebpmf_model/train_2tissues_scores/Brain_Hypothalamus.Brain_Cerebellum-{gene}.K2.ebpmf.rds", gene=['NRXN1', 'NRXN2', 'NRXN3']),
        
rule NRXN:
    input:
        expand("ebpmf_model/train_2tissues/{Tissues}-{gene}.K{K}.ebpmf.rds",
               Tissues = ["Brain_Cortex.Brain_Hippocampus"], # , "Brain_Hypothalamus.Brain_Cerebellum"],
               gene = ['NRXN1', 'NRXN2', 'NRXN3'], K = [2, 3, 5, 10]),

rule collect_counts:
    input:
        expand('Counts/{Tissue}/{gene}.Counts.csv.gz', 
               Tissue = ['Brain_Cortex', 'Brain_Hippocampus', 'Muscle_Skeletal'], gene=genes)
               
rule collect_examples:
    input:
        expand("Counts/{Tissue}/{gene}.Counts.csv.gz", Tissue = ["Brain_Cortex"], gene = ["SRSF3", "RTN2"])
        
rule collect_coverage:
    input:
        expand("coverage/samples/Brain_Cortex/{IndID}.bed.gz", IndID = brain_cortex_samples),
        expand("coverage/samples/Brain_Hippocampus/{IndID}.bed.gz", IndID = brain_hippocampus_samples),
        expand("coverage/samples/Muscle_Skeletal/{IndID}.bed.gz", IndID = muscle_skeletal_samples),
        expand("coverage/samples/Whole_Blood/{IndID}.bed.gz", IndID = whole_blood_samples),
        expand("coverage/samples/Liver/{IndID}.bed.gz", IndID = liver_samples),
        expand("coverage/samples/Lung/{IndID}.bed.gz", IndID = lung_samples),
        expand("coverage/samples/Heart_Atrial_Appendage/{IndID}.bed.gz", IndID = heart_atrial_appendage_samples),
        expand("coverage/samples/Brain_Cerebellum/{IndID}.bed.gz", IndID = brain_cerebellum_samples),
        expand("coverage/samples/Skin_Not_Sun_Exposed_Suprapubic/{IndID}.bed.gz", IndID = skin_not_sun_exposed_suprapubic_samples),
    
rule collect_ebpmf:
    input:
        expand("ebpmf_model/train_5tissues/Brain_Cortex.Brain_Hippocampus.Muscle_Skeletal-{gene}.K{K}.ebpmf.rds",
        gene=genes, K = [2, 3, 5])
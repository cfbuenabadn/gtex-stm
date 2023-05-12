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
        

rule collect_counts:
    input:
        expand('Counts/{Tissue}/{gene}.Counts.csv.gz', 
               Tissue = ['Brain_Cortex', 'Brain_Hippocampus', 'Muscle_Skeletal'], gene=genes)
               
rule collect_examples:
    input:
        expand("Counts/{Tissue}/{gene}.Counts.csv.gz", Tissue = ["Brain_Cortex"], gene = ["SRSF3", "RTN2"])
        
rule Download_By_Parts:
    input:
        expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/Brain_Cortex/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam", 
        IndID = brain_cortex_samples
        ),
        expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/Muscle_Skeletal/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam", 
        IndID = muscle_skeletal_samples
        )
        
rule DifferentialSplicing:
    input:
        "DifferentialSplicing/rMATS/Brain_Cortex_v_Muscle_Skeletal_female_test/output/summary.txt",
        "DifferentialSplicing/rMATS/Brain_Cortex_v_Muscle_Skeletal_train/output/summary.txt",
        "DifferentialSplicing/leafcutter/Brain_Cortex_v_Muscle_Skeletal_female_test/output/leafcutter_cluster_significance.txt",
        "DifferentialSplicing/leafcutter/Brain_Cortex_v_Muscle_Skeletal_train/output/leafcutter_cluster_significance.txt",
        "DifferentialSplicing/leafcutter/juncfiles/Brain_Cortex_v_Muscle_Skeletal/cluster_annotation.tab.gz",
        "DifferentialSplicing/rMATS/Brain_Cortex_v_Muscle_Skeletal_male_test/output/summary.txt",
        expand("DifferentialSplicing/rMATS/negative_tests/{Tissue}_{Sex}_test/output/summary.txt",
            Tissue = ["Brain_Cortex", "Muscle_Skeletal"], Sex=["female", "male"]),
        #"DifferentialSplicing/rMATS/Whole_Blood_v_Liver_female_test/output/summary.txt",
        #"DifferentialSplicing/rMATS/Whole_Blood_v_Liver_train/output/summary.txt",
        #"DifferentialSplicing/leafcutter/Whole_Blood_v_Liver_female_test/output/leafcutter_cluster_significance.txt",
        #"DifferentialSplicing/leafcutter/Whole_Blood_v_Liver_train/output/leafcutter_cluster_significance.txt",
        #"DifferentialSplicing/leafcutter/juncfiles/Whole_Blood_v_Liver/cluster_annotation.tab.gz",
        #"DifferentialSplicing/rMATS/Whole_Blood_v_Muscle_Skeletal_female_test/output/summary.txt",
        #"DifferentialSplicing/rMATS/Whole_Blood_v_Muscle_Skeletal_train/output/summary.txt",
        #"DifferentialSplicing/leafcutter/Whole_Blood_v_Muscle_Skeletal_female_test/output/leafcutter_cluster_significance.txt",
        #"DifferentialSplicing/leafcutter/Whole_Blood_v_Muscle_Skeletal_train/output/leafcutter_cluster_significance.txt",
        #"DifferentialSplicing/leafcutter/juncfiles/Whole_Blood_v_Muscle_Skeletal/cluster_annotation.tab.gz",
        
rule collect_coverage:
    input:
        expand("coverage/bed/Whole_Blood/{IndID}.bed.bgz.tbi", IndID = whole_blood_samples),
        expand("coverage/bed/Liver/{IndID}.bed.bgz.tbi", IndID = liver_samples),
        #expand("coverage/bigwigs/Whole_Blood/{IndID}.bw", IndID = whole_blood_samples),
        #expand("coverage/bigwigs/Liver/{IndID}.bw", IndID = liver_samples),
        expand("/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/Whole_Blood/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc", 
        IndID = whole_blood_samples
        ),
        expand("/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/Liver/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc", 
        IndID = liver_samples
        ),
        #expand("featureCounts/{Tissue}/Counts.txt", Tissue=["Whole_Blood", "Liver"]),
        expand("coverage/bed/Brain_Cortex/{IndID}.bed.bgz.tbi", IndID = brain_cortex_samples),
        expand("coverage/bed/Muscle_Skeletal/{IndID}.bed.bgz.tbi", IndID = muscle_skeletal_samples),
        #expand("coverage/bigwigs/Brain_Cortex/{IndID}.bw", IndID = brain_cortex_samples),
        #expand("coverage/bigwigs/Muscle_Skeletal/{IndID}.bw", IndID = muscle_skeletal_samples),
        expand("/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/Brain_Cortex/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc", 
        IndID = brain_cortex_samples
        ),
        expand("/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/Muscle_Skeletal/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc", 
        IndID = muscle_skeletal_samples
        ),
        #expand("featureCounts/{Tissue}/Counts.txt", Tissue=["Brain_Cortex", "Muscle_Skeletal"])
        
        
        
        
        #expand("coverage/samples/Brain_Hippocampus/{IndID}.bed.gz", IndID = brain_hippocampus_samples),
        #expand("coverage/samples/Muscle_Skeletal/{IndID}.bed.gz", IndID = muscle_skeletal_samples),
        #expand("coverage/samples/Whole_Blood/{IndID}.bed.gz", IndID = whole_blood_samples),
        #expand("coverage/samples/Liver/{IndID}.bed.gz", IndID = liver_samples),
        #expand("coverage/samples/Lung/{IndID}.bed.gz", IndID = lung_samples),
        #expand("coverage/samples/Heart_Atrial_Appendage/{IndID}.bed.gz", IndID = heart_atrial_appendage_samples),
        #expand("coverage/samples/Brain_Cerebellum/{IndID}.bed.gz", IndID = brain_cerebellum_samples),
        #expand("coverage/samples/Skin_Not_Sun_Exposed_Suprapubic/{IndID}.bed.gz", IndID = skin_not_sun_exposed_suprapubic_samples),
        #expand("coverage/samples/Kidney_Cortex/{IndID}.bed.gz", IndID = kidney_cortex_samples),
        
Tissues = ".".join(["Brain_Cerebellum", "Brain_Cortex", "Brain_Hippocampus", 
                            "Heart_Atrial_Appendage", "Kidney_Cortex", "Liver", "Lung", 
                            "Muscle_Skeletal", "Skin_Not_Sun_Exposed_Suprapubic", "Whole_Blood"])
    

rule collect_lm:
    input:
        expand("ebpmf_model/train_2tissues_scores/Brain_Cortex.Muscle_Skeletal/{gene}.K2.ebpmf.rds", gene=genes),
        expand("ebpmf_model/train_2tissues/Brain_Cortex.Muscle_Skeletal/plots/{geneName}.K2.factors.pdf",
        geneName=genes)

rule collect_ebpmf:
    input:
        expand("ebpmf_model/train_10tissues/models/{gene}.K{K}.ebpmf.rds",
        gene=["PKM", "AR"] + quick_test_genes, #genes, 
        K = [2, 3, 4, 5]),
        expand("ebpmf_model/train_10tissues/plots/{geneName}.K{kfactors}.factors.pdf",
        geneName=genes,#["PKM", "AR"] + quick_test_genes, #genes, 
        kfactors = [2, 3, 4, 5]),
        
        

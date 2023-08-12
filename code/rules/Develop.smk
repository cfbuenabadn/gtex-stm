rule DifferentialSplicing:
    input:
        "DifferentialSplicing/rMATS/Brain_Frontal_Cortex_BA9_v_Muscle_Skeletal_female_test/output/summary.txt",
        "DifferentialSplicing/rMATS/Brain_Frontal_Cortex_BA9_v_Muscle_Skeletal_male_test/output/summary.txt",
        "DifferentialSplicing/rMATS/Brain_Frontal_Cortex_BA9_v_Muscle_Skeletal_train/output/summary.txt",
        "DifferentialSplicing/leafcutter/Brain_Frontal_Cortex_BA9_v_Muscle_Skeletal_male_test/output/leafcutter_cluster_significance.txt",
        "DifferentialSplicing/leafcutter/Brain_Frontal_Cortex_BA9_v_Muscle_Skeletal_female_test/output/leafcutter_cluster_significance.txt",
        "DifferentialSplicing/leafcutter/Brain_Frontal_Cortex_BA9_v_Muscle_Skeletal_train/output/leafcutter_cluster_significance.txt",
        "DifferentialSplicing/leafcutter/juncfiles/Brain_Frontal_Cortex_BA9_v_Muscle_Skeletal/cluster_annotation.tab.gz",
        expand("DifferentialSplicing/rMATS/negative_tests/{Tissue}_{Sex}_test/output/summary.txt",
            Tissue = ["Brain_Frontal_Cortex_BA9", "Muscle_Skeletal"], Sex=['female', 'male']),
        #"DifferentialSplicing/rMATS/Cells_EBV-transformed_lymphocytes_v_Cells_Cultured_fibroblasts_female_test/output/summary.txt",
        #"DifferentialSplicing/rMATS/Cells_EBV-transformed_lymphocytes_v_Cells_Cultured_fibroblasts_train/output/summary.txt",
        #"DifferentialSplicing/leafcutter/Cells_EBV-transformed_lymphocytes_v_Cells_Cultured_fibroblasts_female_test/output/leafcutter_cluster_significance.txt",
        #"DifferentialSplicing/leafcutter/Cells_EBV-transformed_lymphocytes_v_Cells_Cultured_fibroblasts_train/output/leafcutter_cluster_significance.txt",
        #"DifferentialSplicing/leafcutter/juncfiles/Cells_EBV-transformed_lymphocytes_v_Cells_Cultured_fibroblasts/cluster_annotation.tab.gz",
        #expand("DifferentialSplicing/rMATS/negative_tests/{Tissue}_female_test/output/summary.txt",
        #    Tissue = ["Cells_EBV-transformed_lymphocytes", "Cells_Cultured_fibroblasts"]),
            
            
rule CollectCoverage:
    input:
        expand("coverage/bed/Brain_Cortex/{IndID}.bed.gz{ext}", IndID = brain_cortex_samples, ext=['', '.tbi']),
        expand("/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/Brain_Cortex/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc", 
        IndID = brain_cortex_samples
        ),
        expand("coverage/bed/Muscle_Skeletal/{IndID}.bed.gz{ext}", IndID = muscle_skeletal_samples, ext=['', '.tbi']),
        expand("/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/Muscle_Skeletal/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc", 
        IndID = muscle_skeletal_samples
        ),
        expand("coverage/bed/Liver/{IndID}.bed.gz{ext}", IndID = liver_samples, ext=['', '.tbi']),
        expand("/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/Liver/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc", 
        IndID = liver_samples
        ),
        expand("coverage/bed/Whole_Blood/{IndID}.bed.gz{ext}", IndID = whole_blood_samples, ext=['', '.tbi']),
        expand("/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/Whole_Blood/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc", 
        IndID = whole_blood_samples
        ),
        expand("coverage/bed/Brain_Hypothalamus/{IndID}.bed.gz{ext}", IndID = brain_hypothalamus_samples, ext=['', '.tbi']),
        expand("/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/Brain_Hypothalamus/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc", 
        IndID = brain_hypothalamus_samples
        ),
        expand("coverage/bed/Brain_Cerebellum/{IndID}.bed.gz{ext}", IndID = brain_cerebellum_samples, ext=['', '.tbi']),
        expand("/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/Brain_Cerebellum/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc", 
        IndID = brain_cerebellum_samples
        ),
        expand("coverage/bed/Brain_Hippocampus/{IndID}.bed.gz{ext}", IndID = brain_hippocampus_samples, ext=['', '.tbi']),
        expand("/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/Brain_Hippocampus/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc", 
        IndID = brain_hippocampus_samples
        ),
        expand("coverage/bed/Cells_Cultured_fibroblasts/{IndID}.bed.gz{ext}", IndID = fibroblast_samples, ext=['', '.tbi']),
        expand("/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/Cells_Cultured_fibroblasts/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc", 
        IndID = fibroblast_samples
        ),
        expand("coverage/bed/Cells_EBV-transformed_lymphocytes/{IndID}.bed.gz{ext}", IndID = LCL_samples, ext=['', '.tbi']),
        expand("/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/Cells_EBV-transformed_lymphocytes/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc", 
        IndID = LCL_samples
        ),
        
rule CollectCoverageByParts:
    input:
        expand("coverage/bed/Brain_Frontal_Cortex_BA9/{IndID}.bed.gz{ext}", IndID = BA9_samples, ext=['', '.tbi']),
        expand("/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/Brain_Frontal_Cortex_BA9/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc", 
        IndID = BA9_samples
        ),
        expand("coverage/bed/Muscle_Skeletal/{IndID}.bed.gz{ext}", IndID = muscle_skeletal_samples, ext=['', '.tbi']),
        expand("/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/Muscle_Skeletal/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc", 
        IndID = muscle_skeletal_samples
        ),
        

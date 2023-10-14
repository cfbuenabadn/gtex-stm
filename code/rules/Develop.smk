rule DifferentialSplicing:
    input:
        "DifferentialSplicing/rMATS/Brain_Frontal_Cortex_BA9_v_Muscle_Skeletal_female_test/output/summary.txt",
        "DifferentialSplicing/rMATS/Brain_Frontal_Cortex_BA9_v_Muscle_Skeletal_male_test/output/summary.txt",
        "DifferentialSplicing/rMATS/Brain_Frontal_Cortex_BA9_v_Muscle_Skeletal_train/output/summary.txt",
        "DifferentialSplicing/leafcutter/Brain_Frontal_Cortex_BA9_v_Muscle_Skeletal_male_test/output/leafcutter_cluster_significance.txt",
        "DifferentialSplicing/leafcutter/Brain_Frontal_Cortex_BA9_v_Muscle_Skeletal_female_test/output/leafcutter_cluster_significance.txt",
        "DifferentialSplicing/leafcutter/Brain_Frontal_Cortex_BA9_v_Muscle_Skeletal_train/output/leafcutter_cluster_significance.txt",
        #"DifferentialSplicing/leafcutter/juncfiles/Brain_Frontal_Cortex_BA9_v_Muscle_Skeletal/cluster_annotation.tab.gz",
        expand("DifferentialSplicing/rMATS/negative_tests/{Tissue}_{Sex}_test/output/summary.txt",
            Tissue = ["Brain_Frontal_Cortex_BA9", "Muscle_Skeletal"], Sex=['female', 'male']),
        "DifferentialSplicing/rMATS/Heart_Atrial_Appendage_v_Lung_female_test/output/summary.txt",
        "DifferentialSplicing/rMATS/Heart_Atrial_Appendage_v_Lung_train/output/summary.txt",
        "DifferentialSplicing/leafcutter/Heart_Atrial_Appendage_v_Lung_female_test/output/leafcutter_cluster_significance.txt",
        "DifferentialSplicing/leafcutter/Heart_Atrial_Appendage_v_Lung_train/output/leafcutter_cluster_significance.txt",
        #"DifferentialSplicing/leafcutter/juncfiles/Heart_Atrial_Appendage_v_Lung/cluster_annotation.tab.gz",
        expand("DifferentialSplicing/rMATS/negative_tests/{Tissue}_female_test/output/summary.txt",
            Tissue = ["Heart_Atrial_Appendage", "Lung"]),
            
            
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
        expand("coverage/bed/Whole_Blood/{IndID}.bed.gz{ext}", IndID = whole_blood_samples, ext=['', '.tbi']),
        expand("/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/Whole_Blood/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc",
        IndID = whole_blood_samples
        ),        
        expand("coverage/bed/Heart_Atrial_Appendage/{IndID}.bed.gz{ext}", IndID = heart_atrial_appendage_samples, ext=['', '.tbi']),
        expand("/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/Heart_Atrial_Appendage/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc",
        IndID = heart_atrial_appendage_samples
        ),
        expand("coverage/bed/Lung/{IndID}.bed.gz{ext}", IndID = lung_samples, ext=['', '.tbi']),
        expand("/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/Lung/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc",
        IndID = lung_samples
        ),
        expand("coverage/bed/Brain_Anterior_cingulate_cortex_BA24/{IndID}.bed.gz{ext}", IndID = BA24_samples, ext=['', '.tbi']),
        expand("/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/Brain_Anterior_cingulate_cortex_BA24/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc", 
        IndID = BA24_samples
        ),
        expand("coverage/bed/Skin_Not_Sun_Exposed_Suprapubic/{IndID}.bed.gz{ext}", IndID = skin_not_sun_exposed_suprapubic_samples, ext=['', '.tbi']),
        expand("/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/Skin_Not_Sun_Exposed_Suprapubic/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc", 
        IndID = skin_not_sun_exposed_suprapubic_samples
        ),
        expand("coverage/bed/Liver/{IndID}.bed.gz{ext}", IndID = liver_samples, ext=['', '.tbi']),
        expand("/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/Liver/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc",
        IndID = liver_samples
        ),
        expand("coverage/bed/Brain_Cortex/{IndID}.bed.gz{ext}", IndID = brain_cortex_samples, ext=['', '.tbi']),
        expand("/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/Brain_Cortex/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc",
        IndID = brain_cortex_samples
        ),
        expand("coverage/bed/Brain_Putamen_basal_ganglia/{IndID}.bed.gz{ext}", IndID = putamen_samples, ext=['', '.tbi']),
        expand("/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/Brain_Putamen_basal_ganglia/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.junc",
        IndID = putamen_samples
        ),

rule CollectCountsCoverage:
    input:
        expand("coverage/counts/{Tissues}/{gene}.csv.gz", Tissues = tissue_list, 
               gene = list(selected_genes.gene))
               
rule run_ebpmf_test:
    input:
        #expand('ebpmf_models/RDS/{gene}.K{k}.ebpmf.rds', gene = list(selected_genes.gene), k = ['3', '5'])
        expand('ebpmf_models/RDS/{gene}.K{K}.ebpmf.rds', gene = (prueba_genes + list(selected_genes.gene)), K=['3', '4', '5'])
        


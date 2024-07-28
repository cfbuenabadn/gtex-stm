zcat gwas/raw_summary_statistics/24076602-GCST005531-EFO_0003885-Build37.f.tsv.gz | tail -n+2 | awk -F'\t' '{print "chr"$1, $2, $2+1, $3, "chr"$1":"$2":"$4":"$5, $6, $7, $8}' OFS='\t' > gwas/temp/GRCh37/multiple_sclerosis.bed

cat gwas/raw_summary_statistics/GCST009325.tsv | tail -n+2 | awk -F'\t' '{print "chr"$1, $2, $2+1, "chr"$1":"$2":"$4":"$3, "chr"$1":"$2":"$4":"$3, $8, $5, $6}' OFS='\t' > gwas/temp/GRCh37/parkinson_disease.bed

zcat gwas/raw_summary_statistics/GCST013196.tsv.gz | tail -n+2 | awk -F'\t' '{print "chr"$1, $2, $2+1, "chr"$1":"$2":"$4":"$3, "chr"$1":"$2":"$4":"$3, $8, $5, $6}' OFS='\t' > gwas/temp/GRCh37/alzheimer_late_onset.bed

zcat gwas/raw_summary_statistics/GCST90011875_buildGRCh37.tsv.gz | tail -n+2 | awk '{print "chr"$2, $3, $3+1, $1, "chr"$2":"$3":"$5":"$4, $10, $7, $8}' OFS='\t' > gwas/temp/GRCh37/cognitive_aspects_of_educational_attainment.bed

zcat gwas/raw_summary_statistics/GCST90012754_buildGRCh37.tsv.gz | tail -n+2 |  awk -F'\t' '{print "chr"$1, $3, $3+1, $2, "chr"$1":"$3":"$4":"$5, $10, $8, $9}' OFS='\t' > gwas/temp/GRCh37/Tau_protein_presence_in_body_fluid.bed

cat gwas/raw_summary_statistics/GCST90014122_buildGRCh37.tsv | tail -n+2 | awk -F'\t' '{print "chr"$1, $3, $3+1, $1, "chr"$2":"$3":"$5":"$4, $8, $6, $7}' OFS='\t' > gwas/temp/GRCh37/lacunar_stroke.bed

zcat gwas/raw_summary_statistics/GCST90018815_buildGRCh37.tsv.gz | tail -n+2 | awk -F'\t' '{print "chr"$1, $2, $2+1, "chr"$1":"$2":"$4":"$3, "chr"$1":"$2":"$4":"$3, $8, $6, $7}' OFS='\t' > gwas/temp/GRCh37/cerebral_aneurysm.bed 

zcat gwas/raw_summary_statistics/GCST90018870_buildGRCh37.tsv.gz | tail -n+2 | awk -F'\t' '{print "chr"$1, $2, $2+1, "chr"$1":"$2":"$4":"$3, "chr"$1":"$2":"$4":"$3, $8, $6, $7}' OFS='\t' > gwas/temp/GRCh37/intracerebral_hemorrhage.bed 

zcat gwas/raw_summary_statistics/GCST90027164_buildGRCh37.tsv.gz | tail -n+2 | awk -F'\t' '{print "chr"$2, $3, $3+1, $1, "chr"$2":"$3":"$5":"$4, $9, $7, $8}' OFS='\t' > gwas/temp/GRCh37/amyotrophic_lateral_sclerosis.bed 

zcat gwas/raw_summary_statistics/GCST90095138_buildGRCh37.tsv.gz | tail -n+2 | awk -F'\t' '{print "chr"$1, $2, $2+1, "chr"$1":"$2":"$4":"$3, "chr"$1":"$2":"$4":"$3, $8, $6, $7}' OFS='\t' > gwas/temp/GRCh37/circulating_levels_of_total-tau.bed 

zcat gwas/raw_summary_statistics/GCST90104539_buildGRCh37.tsv.gz | tail -n+2 | awk -F'\t' '{print "chr"$1, $2, $2+1, "chr"$1":"$2":"$11":"$10, "chr"$1":"$2":"$11":"$10, $6, $4, $5}' OFS='\t' > gwas/temp/GRCh37/stroke.bed 

zcat gwas/raw_summary_statistics/GCST90104540_buildGRCh37.tsv.gz | tail -n+2 | awk -F'\t' '{print "chr"$1, $2, $2+1, "chr"$1":"$2":"$11":"$10, "chr"$1":"$2":"$11":"$10, $6, $4, $5}' OFS='\t' > gwas/temp/GRCh37/ischemic_stroke.bed 
zcat gwas/raw_summary_statistics/GCST90104541_buildGRCh37.tsv.gz | tail -n+2 | awk -F'\t' '{print "chr"$1, $2, $2+1, "chr"$1":"$2":"$11":"$10, "chr"$1":"$2":"$11":"$10, $6, $4, $5}' OFS='\t' > gwas/temp/GRCh37/cardioembolic_stroke.bed 
zcat gwas/raw_summary_statistics/GCST90104542_buildGRCh37.tsv.gz | tail -n+2 | awk -F'\t' '{print "chr"$1, $2, $2+1, "chr"$1":"$2":"$11":"$10, "chr"$1":"$2":"$11":"$10, $6, $4, $5}' OFS='\t' > gwas/temp/GRCh37/large_artery_stroke.bed 
zcat gwas/raw_summary_statistics/GCST90104543_buildGRCh37.tsv.gz | tail -n+2 | awk -F'\t' '{print "chr"$1, $2, $2+1, "chr"$1":"$2":"$11":"$10, "chr"$1":"$2":"$11":"$10, $6, $4, $5}' OFS='\t' > gwas/temp/GRCh37/small_vessel_stroke.bed 
zcat gwas/raw_summary_statistics/GCST90105076_buildGRCh37.tsv.gz | tail -n+2 | awk -F'\t' '{print "chr"$1, $2, $2+1, $9, "chr"$1":"$2":"$3":"$4, $8, $6, $7}' OFS='\t' > gwas/temp/GRCh37/total_cerebral_volume.bed 

tail -n+2 gwas/raw_summary_statistics/GCST90267280_buildGRCh37.tsv | awk -F'\t' '{print "chr"$2, $3, $3+1, $1, "chr"$2":"$3":"$5":"$4, $8, $6, $7}' OFS='\t' > gwas/temp/GRCh37/age_when_finished_full-time_education.bed 

zcat gwas/raw_summary_statistics/Okbay_27225129-EduYears_Main.txt.gz | tail -n+2 | awk -F'\t' '{print "chr"$2, $3, $3+1, $1, "chr"$2":"$3":"$4":"$5, $9, $7, $8}' OFS='\t' > gwas/temp/GRCh37/educational_attainment-years_of_education.bed 
zcat gwas/raw_summary_statistics/bip2021.tsv.gz | tail -n+74 | awk -F'\t' '{print "chr"$1, $2, $2+1, $3, "chr"$1":"$2":"$4":"$5, $8, $6, $7}' OFS='\t' > gwas/temp/GRCh37/bipolar_disorder.bed 
zcat gwas/raw_summary_statistics/panic2019.tsv.gz | tail -n+73 | tail -n+2 | awk -F'\t' '{print "chr"$1, $2, $2+1, $3, "chr"$1":"$2":"$4":"$5, $8, $6, $7}' OFS='\t' > gwas/temp/GRCh37/panic_disorder.bed 
zcat gwas/raw_summary_statistics/ppd2023.tsv.gz | tail -n+2 | awk '{print "chr"$1, $2, $2+1, $3, "chr"$1":"$2":"$4":"$5, $10, $8, $9}' OFS='\t' > gwas/temp/GRCh37/major_depression.bed 
zcat gwas/raw_summary_statistics/scz2022.tsv.gz | tail -n+75 | awk -F'\t' '{print "chr"$1, $3, $3+1, $2, "chr"$1":"$3":"$4":"$5, $11, $9, $10}' OFS='\t' > gwas/temp/GRCh37/schizophrenia.bed 







zcat gwas/raw_summary_statistics/vanderLeeSJ_prePMID_FLV_EAonly.txt.gz | tail -n+2 | awk -F' ' '{print $2, $3, $4, $6, $7, $8, $1}' - | awk -F':' '{print $1, $2}' - | awk -F ' ' '{print "chr"$7, $8, $8+1, $1, "chr"$7":"$8":"$2":"$3, $6, $4, $5}' OFS='\t' > gwas/temp/GRCh37/frontal_lobe_volume.bed 

zcat gwas/raw_summary_statistics/vanderLeeSJ_prePMID_OLV_EAonly.txt.gz | tail -n+2 | awk -F' ' '{print $2, $3, $4, $6, $7, $8, $1}' - | awk -F':' '{print $1, $2}' - | awk -F ' ' '{print "chr"$7, $8, $8+1, $1, "chr"$7":"$8":"$2":"$3, $6, $4, $5}' OFS='\t' > gwas/temp/GRCh37/occipital_lobe_volume.bed 

zcat gwas/raw_summary_statistics/vanderLeeSJ_prePMID_PLV_EAonly.txt.gz | tail -n+2 | awk -F' ' '{print $2, $3, $4, $6, $7, $8, $1}' - | awk -F':' '{print $1, $2}' - | awk -F ' ' '{print "chr"$7, $8, $8+1, $1, "chr"$7":"$8":"$2":"$3, $6, $4, $5}' OFS='\t' > gwas/temp/GRCh37/parietal_lobe_volume.bed 

zcat gwas/raw_summary_statistics/vanderLeeSJ_prePMID_TLV_EAonly.txt.gz | tail -n+2 | awk -F' ' '{print $2, $3, $4, $6, $7, $8, $1}' - | awk -F':' '{print $1, $2}' - | awk -F ' ' '{print "chr"$7, $8, $8+1, $1, "chr"$7":"$8":"$2":"$3, $6, $4, $5}' OFS='\t' > gwas/temp/GRCh37/temporal_lobe_volume.bed 




zcat gwas/raw_summary_statistics/GCST90271608.tsv.gz | tail -n+2 | awk -F'\t' '{print "chr"$1, $2, $2+1, $9, "chr"$1":"$2":"$4":"$3, $8, $5, $6}' OFS='\t' > gwas/temp/GRCh37/epilepsy.bed
zcat gwas/raw_summary_statistics/GCST90271609.tsv.gz | tail -n+2 | awk -F'\t' '{print "chr"$1, $2, $2+1, $9, "chr"$1":"$2":"$4":"$3, $8, $5, $6}' OFS='\t' > gwas/temp/GRCh37/genetic_generalized_epilepsy.bed
zcat gwas/raw_summary_statistics/GCST90271616.tsv.gz | tail -n+2 | awk -F'\t' '{print "chr"$1, $2, $2+1, $9, "chr"$1":"$2":"$4":"$3, $8, $5, $6}' OFS='\t' > gwas/temp/GRCh37/childhood_absence_epilepsy.bed
tail -n+2 gwas/raw_summary_statistics/GCST90000016_buildGRCh37.tsv | awk -F'\t' '{print "chr"$1, $2, $2+1, $3, "chr"$1":"$2":"$5":"$4, $6, $7}' OFS='\t' > gwas/temp/GRCh37/migraine.bed
zcat gwas/raw_summary_statistics/adhd2022.tsv.gz | tail -n+2 | awk '{print "chr"$1, $3, $3+1, $2, "chr"$1":"$3":"$4":"$5, $11, $9, $10}' OFS='\t' > gwas/temp/GRCh37/ADHD.bed
zcat gwas/raw_summary_statistics/asd2019.tsv.gz | tail -n+2 | awk -F'\t' '{print "chr"$1, $3, $3+1, $2, "chr"$1":"$3":"$4":"$5, $9, $7, $8}' OFS='\t' > gwas/temp/GRCh37/autism_spectrum_disorder.bed







zcat gwas/raw_summary_statistics/GCST90027158_buildGRCh38.tsv.gz | tail -n+2 | awk -F'\t' '{print "chr"$3, $4, $4+1, $1, "chr"$3":"$4":"$6":"$5, $2, $11, $12}' OFS='\t' > gwas/temp/GRCh38/alzheimer_disease.bed 
zcat gwas/raw_summary_statistics/GCST90085819_buildGRCh38.tsv.gz | tail -n+2 | awk -F'\t' '{print "chr"$3, $4, $4+1, $1, "chr"$3":"$4":"$6":"$5, $10, $8, $9}' OFS='\t' > gwas/temp/GRCh38/brain_stem_volume.bed 
zcat gwas/raw_summary_statistics/TREM22024.tsv.gz | tail -n+2 | awk -F'\t' '{print "chr"$1, $2, $2+1, $5, "chr"$1":"$2":"$3":"$4, $7, $8, $9}' OFS='\t' > gwas/temp/GRCh38/TREM2_in_cerebrospinal_fluid.bed 
zcat gwas/raw_summary_statistics/GCST90080034_buildGRCh38.tsv.gz | tail -n+2 | awk -F'\t' '{print "chr"$2, $3, $3+1, "chr"$1, "chr"$1, $12, $9}' OFS='\t' > gwas/temp/GRCh38/cerebral_infarction.bed

tail -n+2 gwas/raw_summary_statistics/GCST90129600_buildGRCh38.tsv | awk -F'\t' '{print "chr"$1, $2, $2+1, "chr"$1":"$2":"$4":"$3, "chr"$1":"$2":"$4":"$3, $7, $6, $8, $5}' OFS='\t' > gwas/temp/GRCh38/cerebrospinal_fluid_p-tau_levels.bed 

zcat gwas/raw_summary_statistics/GCST90134630_buildGRCh38.tsv.gz | tail -n+2 | awk -F'\t' '{print "chr"$1, $2, $2+1, "chr"$1":"$2":"$4":"$3, "chr"$1":"$2":"$4":"$3, $7, $6, $8, $5}' OFS='\t' > gwas/temp/GRCh38/cerebrospinal_fluid_p-tau_levels_in_abnormal_amyloid_levels.bed 

zcat gwas/raw_summary_statistics/GCST90134631_buildGRCh38.tsv.gz | tail -n+2 | awk -F'\t' '{print "chr"$1, $2, $2+1, "chr"$1":"$2":"$4":"$3, "chr"$1":"$2":"$4":"$3, $7, $6, $8, $5}' OFS='\t' > gwas/temp/GRCh38/cerebrospinal_fluid_p-tau_levels_in_APOE_e4_non-carriers.bed 

zcat gwas/raw_summary_statistics/GCST90134632_buildGRCh38.tsv.gz | tail -n+2 | awk -F'\t' '{print "chr"$1, $2, $2+1, "chr"$1":"$2":"$4":"$3, "chr"$1":"$2":"$4":"$3, $7, $6, $8, $5}' OFS='\t' > gwas/temp/GRCh38/cerebrospinal_fluid_p-tau_levels_in_APOE_e4_carriers.bed 

zcat gwas/raw_summary_statistics/GCST90134633_buildGRCh38.tsv.gz | tail -n+2 | awk -F'\t' '{print "chr"$1, $2, $2+1, "chr"$1":"$2":"$4":"$3, "chr"$1":"$2":"$4":"$3, $7, $6, $8, $5}' OFS='\t' > gwas/temp/GRCh38/cerebrospinal_fluid_p-tau_levels_in_normal_amyloid_levels.bed 






zcat gwas/raw_summary_statistics/AD_sumstats_Jansenetal_2019sept.txt.gz | tail -n+2 | awk '{print "chr"$2, $3, $3+1, $6, "chr"$2":"$3":"$4":"$5, $8, $13, $14}' OFS='\t' > gwas/temp/GRCh37/alzheimer_disease_2.bed 

tail -n+2 gwas/raw_summary_statistics/4_UK_Biobank_IGAP_17May2018 | awk '{print "chr"$8, $9, $9+1, "chr"$8":"$9":"$2":"$3, "chr"$8":"$9":"$2":"$3, $6, $4, $5}' OFS='\t' > gwas/temp/GRCh37/alzheimer_disease_3.bed 

tail -n+2 gwas/raw_summary_statistics/GCST011365_buildGRCh37.tsv | awk '{print "chr"$3, $4, $4+1, $1, "chr"$3":"$4":"$6":"$5, $2, $11, $12}' OFS='\t' > gwas/temp/GRCh37/myocardial_infarction.bed

zcat gwas/raw_summary_statistics/29892015-GCST006061-EFO_0000275-build37.f.tsv.gz | tail -n+2 | awk '{print "chr"$4, $5, $5+1, $1, "chr"$4":"$5":"$3":"$2, $8, $6, $7}' OFS='\t' > gwas/temp/GRCh37/atrial_fibrillation.bed

tail -n+2 gwas/raw_summary_statistics/GCST90132314_buildGRCh37.tsv | awk '{print "chr"$2, $3, $3+1, "chr"$2":"$3":"$5":"$4, "chr"$2":"$3":"$5":"$4, $1, $8, $9}' OFS='\t' > gwas/temp/GRCh37/coronary_artery_disease.bed

zcat gwas/raw_summary_statistics/GCST90162626_buildGRCh37.tsv.gz | tail -n+2 | awk '{print "chr"$2, $3, $3+1, $1, "chr"$2":"$3":"$5":"$4, $9, $7, $8}' OFS='\t' > gwas/temp/GRCh37/heart_failure.bed

zcat gwas/raw_summary_statistics/GCST90267355.tsv.gz | tail -n+2 | awk '{print "chr"$3, $4, $4+1, $1, "chr"$3":"$4":"$6":"$5, $2, $8, $9}' OFS='\t' > gwas/temp/GRCh37/posterior_thigh_muscle_fat_infiltration_percentage.bed

zcat gwas/raw_summary_statistics/GCST90137411_buildGRCh37.tsv.gz | tail -n+2 | awk '{print "chr"$2, $3, $3+1, $1, "chr"$2":"$3":"$5":"$4, $12, $10, $11}' OFS='\t' > gwas/temp/GRCh37/basal_cell_carcinoma.bed

zcat gwas/raw_summary_statistics/34594039-GCST90018921-EFO_0004198-Build37.f.tsv.gz | tail -n+2 | awk '{print "chr"$1, $2, "chr"$1":"$2":"$4":"$3, "chr"$1":"$2":"$4":"$3, $8, $6, $7}' OFS='\t' > gwas/temp/GRCh37/skin_cancer.bed




 


 




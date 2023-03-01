rule gather_counts:
    input:
        expand("Counts/{Tissue}/{gene}.Counts.csv.gz", 
               Tissue = ["Brain_Cortex", "Muscle_Skeletal", "Liver", "Whole_Blood"],
               gene=genes)

rule collect_coloc:
    input:
        expand("hyprcoloc/Results/snmf_{K}.ForGWASColoc/GWAS_GTEx_01/results.txt.gz", K = [3, 5, 10]),
        "hyprcoloc/Results/transcripts_10.ForGWASColoc/GWAS_GTEx_01/results.txt.gz"


rule CollectPi1Tables:
    input:
        expand("pi1/pi1_tables/Pi1.snmf_{K}..txt", #SexGroup = ["", ".male", ".female"],
               K = [3, 5, 10]),
        "pi1/pi1_tables/Pi1.transcripts_10..txt"
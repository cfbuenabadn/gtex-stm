
# use tabix genewise qtltools output to ascertain p value in trait2 for top SNP determined from trait1 discoveries.
# There are so many potential trait pairs, so i am going to break this down into chunks to parralelize
rule CalculatePi1_GetTraitPairs_AllTraits:
    input:
        Nominal = expand("QTLs/{tissue_id}/multi_tissue{{SexGroup}}.snmf_{{K}}.NominalPass.txt.tabix.gz", tissue_id = tissue_sub_list),
        Permutation = expand("QTLs/{tissue_id}/multi_tissue{{SexGroup}}.snmf_{{K}}.PermutationPass.txt.gz", tissue_id = tissue_sub_list)
    params:
        NumChunks = NumPvalsForPi1Chunks
    output:
        expand("pi1/PairwiseTraitsToCompare/Pairwise{{SexGroup}}.snmf_{{K}}.{n}.txt.gz", n=range(1, 1+NumPvalsForPi1Chunks))
    conda:
        "../envs/r_essentials.yml"
    wildcard_constraints:
        SexGroup = "|.female|.male",
        K = "2|3|4|5|10"
    resources:
        mem_mb = 32000
    log:
        "logs/CalculatePi1_GetTraitPairs_AllTraits.{SexGroup}.snmf_{K}.log"
    shell:
        """
        Rscript scripts/CalculatePi1_GetTraitPairs_AllTraits.R {params.NumChunks} pi1/PairwiseTraitsToCompare/Pairwise{wildcards.SexGroup}.snmf_{wildcards.K}. {input.Permutation} &> {log}
        """

rule GatherChunks:
    input:
        expand("pi1/PairwiseTraitsToCompare/Pairwise{SexGroup}.snmf_{K}.{n}.txt.gz", n=range(1, 1+NumPvalsForPi1Chunks), SexGroup = ["", ".male", ".female"],
               K = [2, 3, 4, 5, 10])

rule GetPvalsForPi1AllTraitPairs:
    input:
        TraitsToCompare = "pi1/PairwiseTraitsToCompare/Pairwise{SexGroup}.snmf_{K}.{chunk}.txt.gz",
        tabix_QTLsOut = expand("QTLs/{tissue_id}/multi_tissue{{SexGroup}}.snmf_{{K}}.NominalPass.txt.tabix.gz", tissue_id=tissue_sub_list)
    output:
        "pi1/PairwiseTraitsToCompare/P{SexGroup}.snmf_{K}.{chunk}.txt.gz"
    log:
        "logs/GetPvalsForPi1AllTraitPairs/Ascertainment{SexGroup}.snmf_{K}.{chunk}.log"
    wildcard_constraints:
        SexGroup = "|.female|.male",
        K = "2|3|4|5|10"
    shell:
        """
        python scripts/CalculatePi1_GetAscertainmentP_AllPairs.py {input.TraitsToCompare} {output} {input.tabix_QTLsOut} &> {log}
        """

rule GatherPvalsForPi1AllTraitPairs:
    input:
        expand("pi1/PairwiseTraitsToCompare/P{{SexGroup}}.snmf_{{K}}.{chunk}.txt.gz", chunk=range(1, NumPvalsForPi1Chunks+1))
    output:
        "pi1/input_tables/PairwisePi1Traits{SexGroup}.snmf_{K}.P.all.txt.gz"
    wildcard_constraints:
        SexGroup = "|.female|.male",
        K = "2|3|4|5|10"
    shell:
        """
        cat <(zcat {input[0]} | head -1) <(zcat {input} | grep -v -P '^PC1\\t') | gzip - > {output}
        """

rule CalculatePi1Table:
    input:
        "pi1/input_tables/PairwisePi1Traits{SexGroup}.snmf_{K}.P.all.txt.gz"
    output:
        "pi1/pi1_tables/Pi1{SexGroup}.snmf_{K}.txt"
    log:
        "logs/pi1_tables/Pi1{SexGroup}.snmf_{K}.log"
    wildcard_constraints:
        SexGroup = "|.female|.male",
        K = "2|3|4|5|10"
    shell:
        """
        Rscript scripts/CalculatePi1Statistics.R {input} {output} &> {log}
        """

rule CollectPi1Tables:
    input:
        expand("pi1/pi1_tables/Pi1.snmf_{K}.txt", #SexGroup = ["", ".male", ".female"],
               K = [3, 5, 10])
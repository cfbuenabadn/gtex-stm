# use tabix genewise qtltools output to ascertain p value in trait2 for top SNP determined from trait1 discoveries.
# There are so many potential trait pairs, so i am going to break this down into chunks to parralelize

rule CalculatePi1_GetTraitPairs_AllTraits:
    input:
        Nominal = expand("QTLs/GTEx_10/{tissue_id}/{{annotation}}_{{K}}{{SexGroup}}{{QTLsRun}}NominalPass.txt.tabix.gz", tissue_id = tissue_sub_list),
        Permutation = expand("QTLs/GTEx_10/{tissue_id}/{{annotation}}_{{K}}{{SexGroup}}{{QTLsRun}}PermutationPass.txt.gz", tissue_id = tissue_sub_list)
    params:
        NumChunks = NumPvalsForPi1Chunks
    output:
        expand("pi1/PairwiseTraitsToCompare/Pairwise{{SexGroup}}.{{annotation}}_{{K}}{{QTLsRun}}{n}.txt.gz", n=range(1, 1+NumPvalsForPi1Chunks))
    conda:
        "../envs/r_essentials.yml"
    wildcard_constraints:
        SexGroup = "|.female|.male",
        K = "2|3|4|5|10",
        annotation = 'transcripts|snmf',
        QTLsRun = '.Grouped|.'
    resources:
        mem_mb = 32000
    log:
        "logs/CalculatePi1_GetTraitPairs_AllTraits.{SexGroup}.{annotation}_{K}{QTLsRun}.log"
    shell:
        """
        Rscript scripts/CalculatePi1_GetTraitPairs_AllTraits.R {params.NumChunks} pi1/PairwiseTraitsToCompare/Pairwise{wildcards.SexGroup}.{wildcards.annotation}_{wildcards.K}. {input.Permutation} &> {log}
        """

#rule GatherChunks:
#    input:
#        expand("pi1/PairwiseTraitsToCompare/Pairwise.snmf_{K}.{n}.txt.gz", n=range(1, 1+NumPvalsForPi1Chunks),
#        expand("pi1/PairwiseTraitsToCompare/Pairwise.transcripts_10.{n}.txt.gz", n=range(1, 1+NumPvalsForPi1Chunks))



rule GetPvalsForPi1AllTraitPairs:
    input:
        TraitsToCompare = "pi1/PairwiseTraitsToCompare/Pairwise{SexGroup}.{annotation}_{K}{QTLsRun}{chunk}.txt.gz",
        tabix_QTLsOut = expand("QTLs/GTEx_10/{tissue_id}/{{annotation}}_{{K}}{{SexGroup}}{{QTLsRun}}NominalPass.txt.tabix.gz", tissue_id=tissue_sub_list)
    output:
        "pi1/PairwiseTraitsToCompare/P{SexGroup}.{annotation}_{K}{QTLsRun}{chunk}.txt.gz"
    log:
        "logs/GetPvalsForPi1AllTraitPairs/Ascertainment{SexGroup}.{annotation}_{K}{QTLsRun}.{chunk}.log"
    wildcard_constraints:
        SexGroup = "|.female|.male",
        K = "2|3|4|5|10",
        annotation = 'transcripts|snmf',
        QTLsRun = ".Grouped|."
    resources:
        mem_mb = 32000
    shell:
        """
        python scripts/CalculatePi1_GetAscertainmentP_AllPairs.py {input.TraitsToCompare} {output} {input.tabix_QTLsOut} &> {log}
        """

rule GatherPvalsForPi1AllTraitPairs:
    input:
        expand("pi1/PairwiseTraitsToCompare/P{{SexGroup}}.{{annotation}}_{{K}}{{QTLsRun}}{chunk}.txt.gz", chunk=range(1, NumPvalsForPi1Chunks+1))
    output:
        "pi1/input_tables/PairwisePi1Traits{SexGroup}.{annotation}_{K}{QTLsRun}P.all.txt.gz"
    wildcard_constraints:
        SexGroup = "|.female|.male",
        K = "2|3|4|5|10",
        annotation = 'transcripts|snmf',
        QTLsRun = ".Grouped|."
    shell:
        """
        cat <(zcat {input[0]} | head -1) <(zcat {input} | grep -v -P '^PC1\\t') | gzip - > {output}
        """

rule CalculatePi1Table:
    input:
        "pi1/input_tables/PairwisePi1Traits{SexGroup}.{annotation}_{K}{QTLsRun}P.all.txt.gz"
    output:
        "pi1/pi1_tables/Pi1{SexGroup}.{annotation}_{K}{QTLsRun}txt"
    log:
        "logs/pi1_tables/Pi1{SexGroup}.{annotation}_{K}{QTLsRun}.log"
    wildcard_constraints:
        SexGroup = "|.female|.male",
        K = "2|3|4|5|10",
        annotation = "transcripts|snmf",
        QTLsRun = ".Grouped|."
    resources:
        mem_mb = 32000
    shell:
        """
        Rscript scripts/CalculatePi1Statistics.R {input} {output} &> {log}
        """

               
              
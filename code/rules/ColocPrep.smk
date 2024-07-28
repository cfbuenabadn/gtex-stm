rule SplitAndCombineSummaryStatsPerLocus:
    """
    since hyprcoloc will require reading in summary stat data in R, it will be
    more convenient if the summary stats were organized into smaller files so
    we don't have to read in huge files into R. Also, it is sort of pointless
    to attempt colocalization on molecular traits that don't have an QTL
    signal. Therefore, this script reads the QTLtools output files for multiple
    phenotypes, and writes out the necessary summary stats for all phenotypes
    to new smaller files (one file for each gwas trait or gene) if they pass
    some minimum Pvalue treshold from the QTLtools permutation test.
    """
    input:
        zip( expand(
            "QTLs/GTEx_10/{tissue_id}/{{annotation}}_{{K}}{{SexGroup}}.ForGWASColoc.PermutationPass.txt.gz",
            tissue_id=tissue_sub_list,
        ),
        expand(
            "QTLs/GTEx_10/{tissue_id}/{{annotation}}_{{K}}{{SexGroup}}.ForGWASColoc.NominalPass.Chunks/{{n}}.txt",
            tissue_id=tissue_sub_list,
        )) 
    output:
        directory("hyprcoloc/{annotation}_{K}{SexGroup}.LociWiseSummaryStatsInput/Chunks/ForGWASColoc/{n}")
    log:
        "logs/SplitAndCombineSummaryStatsPerGene.{annotation}_{K}{SexGroup}.ForGWASColoc.{n}.log"
    resources:
        mem_mb = much_more_mem_after_first_attempt
    params:
        MinNominalP = 0.01
    wildcard_constraints:
        SexGroup = "|.female|.male",
        K = '2|3|4|5|10',
        annotation = 'snmf|transcripts'
    shell:
        """
        python scripts/CombineAndSplitSummaryStatsForGWASColoc.py {output}/ {params.MinNominalP} {input} &> {log}
        """

rule GatherSummaryStatsFileChunks:
    input:
        expand("hyprcoloc/{{annotation}}_{{K}}{{SexGroup}}.LociWiseSummaryStatsInput/Chunks/ForGWASColoc/{n}", n=range(1, 1+N_PermutationChunks) )
    output:
        directory("hyprcoloc/{annotation}_{K}{SexGroup}.LociWiseSummaryStatsInput/ForGWASColoc")
    log:
        "logs/GatherSummaryStatsFileChunks/{annotation}_{K}{SexGroup}.log"
    resources:
        mem_mb = much_more_mem_after_first_attempt
    wildcard_constraints:
        SexGroup = "|.female|.male",
        K = '2|3|4|5|10',
        annotation = 'snmf|transcripts'
    shell:
        """
        python scripts/MergeSummaryStatChunks.py {output} {input} &> {log}
        """

rule GatherHyprcolocInput:
    input:
        "hyprcoloc/transcripts_10.LociWiseSummaryStatsInput/ForGWASColoc",
        #"hyprcoloc/snmf_5.LociWiseSummaryStatsInput/ForGWASColoc",
        #"hyprcoloc/snmf_10.LociWiseSummaryStatsInput/ForGWASColoc",



rule create_gwascoloc_bash_scripts:
    output:
        expand("hyprcoloc/Results/{{annotation}}_{{K}}.ForGWASColoc/{{ColocName}}/Chunks/{n}.sh", n=range(0,config["gwas_coloc_chunks"]))
    wildcard_constraints:
        FeatureCoordinatesRedefinedFor = "ForGWASColoc",
        ColocName = '|'.join(colocs_gwas.index)
    params:
        PhenotypesToColoc = GetMolPhenotypesToColoc
    run:
        from itertools import cycle
        gwas_bashscript_pairs = zip(gwas_traits_for_coloc, cycle(output))
        for accession, out_f in gwas_bashscript_pairs:
            with open(out_f, 'a') as f:
                _ = f.write(f'Rscript scripts/hyprcoloc_gwas2.R hyprcoloc/{wildcards.annotation}_{wildcards.K}.LociWiseSummaryStatsInput/ForGWASColoc/{accession}.txt.gz gwas_summary_stats/StatsForColoc/{accession}.standardized.txt.gz {out_f.rstrip(".sh")}.txt.gz "{params.PhenotypesToColoc}"\n')
        # If there are more output files than accession numbers, the extra
        # output files won't get made in the previous loop and snakemake will
        # complain of missing output files. as a fail safe, let's append to
        # each file in output, in effect making an empty file if a file wasn't
        # made in the for loop above
        for f in output:
            open(f, 'a').close()

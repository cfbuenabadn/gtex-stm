rule ComputeMatrixForCoverageMetaplots:
    input:
        bigwigs = "CoveragePlots/bigwigs/{Tissue}.bw",
        bed = expand("CoveragePlots/bed_files/introns/{ri_group}.retained_introns.bed.gz", ri_group = ['snmf_and_gencode', 'snmf_only', 'gencode_only']) 
    output:
        mat = "CoveragePlots/matrices/{Tissue}.mat.gz"
    wildcard_constraints:
        Tissue = "|".join(tissue_sub_list)
    log:
        "logs/Metaplots/ComputeMatrix.{Tissue}.log"
    resources:
        mem_mb = 12000
    shell:
        """
        computeMatrix scale-regions -S {input.bigwigs} -R {input.bed} -b 100 -a 100 --startLabel "5' splice site" --endLabel "3' splice site" -o {output} &> {log}
        """
        

rule plotHeatmapCoverageMetaplots:
    input:
        "CoveragePlots/matrices/{Tissue}.mat.gz"
    output:
        "CoveragePlots/plots/{Tissue}.png"
    wildcard_constraints:
        Tissue = "|".join(tissue_sub_list)
    log:
        "logs/Metaplots/PlotHeatmap.{Tissue}.log"
    shell:
        """
        plotHeatmap -m {input} -o {output} --averageTypeSummaryPlot median --regionsLabel "sNMF and Gencode" "sNMF only" "Gencode only" --heatmapHeight 10 &>> {log}
        """

rule collect_heatmaps:
    input:
        expand('CoveragePlots/plots/{Tissue}.png', Tissue=tissue_sub_list)




rule DownloadFromGTEx_Subset:
    input:
        manifest = 'gtex-download/subset/files/tissue-manifest.json'
    output:
        expand(
        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/subset/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam", 
        IndID = subsamples
        )
    log:
        'logs/download_subsamples.log' 
    resources:
        mem_mb = 42000
    shell:
        """
        (./gen3-client download-multiple --profile=AnVIL --manifest={input.manifest} --download-path=/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/subset/bams/ --protocol=s3) &> {log}
        """

#rule IndexBam_BW:
#    input:
#        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/subset/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam"
#    output:
#        "/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/subset/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.bai"
#    log:
#        "logs/indexbam/subset.{IndID}.log"
#    resources:
#        mem_mb = 24000
##    shell:
#        """
#        samtools index {input} {output} > {log}
#        """

def get_bam_for_tissue_merge(wildcards):
    tissue_sample_list = list(gtex_samples.loc[(gtex_samples.tissue_id == wildcards.Tissue) & (gtex_samples.group == 'train')].index[[0,25, 50, 75, 99]])
    bam_files = expand("/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/subset/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam",
                        IndID = tissue_sample_list)
    return bam_files

def get_bai_for_tissue_merge(wildcards):
    tissue_sample_list = list(gtex_samples.loc[(gtex_samples.tissue_id == wildcards.Tissue) & (gtex_samples.group == 'train')].index[[0,25, 50, 75, 99]])
    bai_files = expand("/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex-download/subset/bams/{IndID}.Aligned.sortedByCoord.out.patched.md.bam.bai",
                        IndID = tissue_sample_list)
    return bai_files

rule MergeBamsForBigWig:
    input:
        bam = get_bam_for_tissue_merge,
        bai = get_bai_for_tissue_merge
    output:
        bam = temp("CoveragePlots/merged_bams/{Tissue}.bam"),
        bam_sorted = "CoveragePlots/merged_bams/{Tissue}.sorted.bam",
        bai = "CoveragePlots/merged_bams/{Tissue}.sorted.bam.bai"
    log:
        "/scratch/midway3/cnajar/logs/mergebams/{Tissue}.log"
    resources:
        mem_mb = 24000
    wildcard_constraints:
        Tissue = '|'.join(tissue_sub_list)
    shell:
        """
        samtools merge -o {output.bam} {input.bam} &> {log};
        samtools sort -o {output.bam_sorted} {output.bam} &>> {log};
        (samtools index {output.bam_sorted} {output.bai} > {log}) &>> {log}
        """

rule BamToBigwig:
    input:
        bam = "CoveragePlots/merged_bams/{Tissue}.sorted.bam",
        bai = "CoveragePlots/merged_bams/{Tissue}.sorted.bam.bai",
    output:
        "CoveragePlots/bigwigs/{Tissue}.bw"
    log:
        "/scratch/midway3/cnajar/logs/bigwigs/{Tissue}.log"
    resources:
        mem_mb = 24000
    wildcard_constraints:
        Tissue = '|'.join(tissue_sub_list)
    shell:
        """
        bamCoverage --binSize 1 --normalizeUsing CPM --effectiveGenomeSize 2913022398 -b {input.bam} -o {output} &> {log}
        """

rule collect_bw:
    input:
        expand('CoveragePlots/bigwigs/{Tissue}.bw', Tissue = tissue_sub_list)


rule GetIntronBeds:
    input:
        'Annotations/gencode.v44.primary_assembly.exons.sorted.bed.gz',
        'ebpmf_models/filtered/snmf_10/tables/snmf.merged_isoforms.exons.sorted.bed.gz',
        'ebpmf_models/filtered/snmf_10/tables/annotated.snmf.merged_isoforms.tab.gz'
    output:
        'CoveragePlots/bed_files/introns/snmf.retained_introns.bed.gz',
        'CoveragePlots/bed_files/introns/gencode.retained_introns.bed.gz',
        'CoveragePlots/bed_files/introns/gencode.appris_introns.bed.gz'
    log:
        "/scratch/midway3/cnajar/logs/bigwigs/get_introns_bed.log"
    resources:
        mem_mb = 12000
    shell:
        """
        python scripts/get_retained_introns.py &> {log}
        """

rule GetInputForIntronMetaplot:
    input:
        snmf = 'CoveragePlots/bed_files/introns/snmf.retained_introns.bed.gz',
        gencode = 'CoveragePlots/bed_files/introns/gencode.retained_introns.bed.gz',
        gencode_exons = 'Annotations/gencode.v44.primary_assembly.exons.sorted.bed.gz',
        snmf_exons = 'ebpmf_models/filtered/snmf_10/tables/snmf.merged_isoforms.exons.sorted.bed.gz',
        appris = 'CoveragePlots/bed_files/introns/gencode.appris_introns.bed.gz'
    output:
        snmf_only_temp = temp('CoveragePlots/bed_files/introns/snmf_only_temp.retained_introns.bed.gz'),
        snmf_only = 'CoveragePlots/bed_files/introns/snmf_only.retained_introns.bed.gz',
        gencode_only_temp = temp('CoveragePlots/bed_files/introns/gencode_only_temp.retained_introns.bed.gz'),
        gencode_only = 'CoveragePlots/bed_files/introns/gencode_only.retained_introns.bed.gz',
        snmf_gencode = 'CoveragePlots/bed_files/introns/snmf_and_gencode.retained_introns.bed.gz',
        appris = 'CoveragePlots/bed_files/introns/appris_introns.bed.gz',
    log:
        "/scratch/midway3/cnajar/logs/bigwigs/get_introns_bed_for_metaplots.log"
    resources:
        mem_mb = 12000
    shell:
        """
        (bedtools intersect -v -a {input.snmf} -b {input.gencode} | bedtools sort -i - > {output.snmf_only_temp}) &>> {log};
        (bedtools intersect -v -f 0.01 -a {output.snmf_only_temp} -b {input.gencode_exons} > {output.snmf_only}) &>> {log};
        (bedtools intersect -v -a {input.gencode} -b {input.snmf} | bedtools sort -i - > {output.gencode_only_temp}) &>> {log};
        (bedtools intersect -v -f 0.01 -a {output.gencode_only_temp} -b {input.snmf_exons} > {output.gencode_only}) &>> {log};
        (bedtools intersect -s -F 1 -f 1 -wa -a {input.snmf} -b {input.gencode} | bedtools sort -i - > {output.snmf_gencode}) &>> {log};
        (bedtools intersect -v -f 0.01 -F 0.01 -a {input.appris} -b {input.snmf_exons} | bedtools intersect -v -f 0.01 -F 0.01 -a - -b {input.gencode_exons} | bedtools sort -i - > {output.appris}) &>> {log}
        """


rule GatherIntronMetaplotInput:
    input:
        'CoveragePlots/bed_files/introns/snmf_only.retained_introns.bed.gz',
        'CoveragePlots/bed_files/introns/gencode_only.retained_introns.bed.gz',
        'CoveragePlots/bed_files/introns/snmf_and_gencode.retained_introns.bed.gz',
        'CoveragePlots/bed_files/introns/appris_introns.bed.gz'


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




rule ComputeMatrixForCoverageMetaplots_sep:
    input:
        bigwigs = "CoveragePlots/bigwigs/{Tissue}.bw",
        bed = "CoveragePlots/bed_files/introns/{ri_group}.bed.gz" 
    output:
        mat = "CoveragePlots/matrices/{Tissue}.{ri_group}.mat.gz"
    wildcard_constraints:
        Tissue = "|".join(tissue_sub_list),
        ri_group = '|'.join(['snmf_and_gencode.retained_introns', 'snmf_only.retained_introns', 'gencode_only.retained_introns', 'appris_introns'])
    log:
        "logs/Metaplots/ComputeMatrix.{Tissue}.{ri_group}.log"
    resources:
        mem_mb = 12000
    shell:
        """
        computeMatrix scale-regions -S {input.bigwigs} -R {input.bed} -b 50 -a 50 --startLabel "5' splice site" --endLabel "3' splice site" -o {output} &> {log}
        """


rule collect_metaplot_matrices:
    input:
        expand('CoveragePlots/matrices/{Tissue}.{ri_group}.mat.gz', Tissue=tissue_sub_list, ri_group = ['snmf_and_gencode.retained_introns', 'snmf_only.retained_introns', 'gencode_only.retained_introns', 'appris_introns'])



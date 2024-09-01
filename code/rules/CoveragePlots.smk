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
        'CoveragePlots/bed_files/introns/vastdb.retained_introns.bed.gz',
        'CoveragePlots/bed_files/introns/gencode.appris_introns.bed.gz',
        'CoveragePlots/bed_files/introns/snmf.flanking_exons_1.bed',
        'CoveragePlots/bed_files/introns/gencode.flanking_exons_1.bed',
        'CoveragePlots/bed_files/introns/vastdb.flanking_exons_1.bed',
        'CoveragePlots/bed_files/introns/appris_flanking_exons_1.bed',
        'CoveragePlots/bed_files/introns/snmf.flanking_exons_2.bed',
        'CoveragePlots/bed_files/introns/gencode.flanking_exons_2.bed',
        'CoveragePlots/bed_files/introns/vastdb.flanking_exons_2.bed',
        'CoveragePlots/bed_files/introns/appris_flanking_exons_2.bed'
    log:
        "/scratch/midway3/cnajar/logs/bigwigs/get_introns_bed.log"
    resources:
        mem_mb = 24000
    shell:
        """
        python scripts/get_retained_introns.py &> {log}
        """

rule GetInputForIntronMetaplot:
    input:
        snmf = 'CoveragePlots/bed_files/introns/snmf.retained_introns.bed.gz',
        gencode = 'CoveragePlots/bed_files/introns/gencode.retained_introns.bed.gz',
        gencode_exons = 'Annotations/gencode.v44.primary_assembly.exons.bed.gz', #'Annotations/gencode.v44.primary_assembly.exons.sorted.bed.gz',
        snmf_exons = 'ebpmf_models/filtered/snmf_10/tables/snmf.merged_isoforms.exons.sorted.bed.gz',
        appris = 'CoveragePlots/bed_files/introns/gencode.appris_introns.bed.gz'
    output:
        snmf_only_temp = temp('CoveragePlots/bed_files/introns/snmf_only_temp.retained_introns.bed'),
        snmf_only = 'CoveragePlots/bed_files/introns/snmf_only.retained_introns.bed',
        gencode_only_temp = temp('CoveragePlots/bed_files/introns/gencode_only_temp.retained_introns.bed'),
        gencode_only = 'CoveragePlots/bed_files/introns/gencode_only.retained_introns.bed',
        snmf_gencode = 'CoveragePlots/bed_files/introns/snmf_and_gencode.retained_introns.bed',
        appris = 'CoveragePlots/bed_files/introns/appris_introns.bed'
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

rule get_bed_for_all_snmf_metaplot:
    input:
        'CoveragePlots/bed_files/introns/snmf.retained_introns.bed.gz'
    output:
        'CoveragePlots/bed_files/introns/snmf.retained_introns.bed'
    log:
        "/scratch/midway3/cnajar/logs/get_unzipped_bed_snmf_introns.log"
    resources:
        mem_mb = 12000
    shell:
        """
        (zcat {input} > {output}) &> {log}
        """

rule GatherIntronMetaplotInput:
    input:
        'CoveragePlots/bed_files/introns/snmf_only.retained_introns.bed',
        'CoveragePlots/bed_files/introns/gencode_only.retained_introns.bed',
        'CoveragePlots/bed_files/introns/snmf_and_gencode.retained_introns.bed',
        'CoveragePlots/bed_files/introns/appris_introns.bed'


def getComputeMatrixParams(wildcards):
    if 'introns' in wildcards.ri_group:
        return '--regionBodyLength 1000'
    else:
        return '--regionBodyLength 200'


rule ComputeMatrixForCoverageMetaplots_sep:
    input:
        bigwigs = "CoveragePlots/bigwigs/{Tissue}.bw",
        bed = "CoveragePlots/bed_files/introns/{ri_group}.bed" 
    output:
        mat = "CoveragePlots/matrices/{Tissue}.{ri_group}.mat.gz"
    wildcard_constraints:
        Tissue = "|".join(tissue_sub_list),
        ri_group = '|'.join(['snmf_and_gencode.retained_introns', 'snmf_only.retained_introns', 'gencode_only.retained_introns', 'appris_introns',
                             'snmf.flanking_exons_1', 'gencode.flanking_exons_1', 'appris_flanking_exons_1',
                             'snmf.flanking_exons_2', 'gencode.flanking_exons_2', 'appris_flanking_exons_2',
                             'vastdb.flanking_exons_1', 'vastdb.flanking_exons_2', 'snmf.retained_introns'])
    log:
        "logs/Metaplots/ComputeMatrix.{Tissue}.{ri_group}.log"
    resources:
        mem_mb = 12000
    params:
        getComputeMatrixParams
    shell:
        """
        computeMatrix scale-regions -S {input.bigwigs} -R {input.bed} {params} -b 50 -a 50 --startLabel "5' splice site" --endLabel "3' splice site" -o {output} &> {log}
        """


rule collect_metaplot_matrices:
    input:
        expand('CoveragePlots/matrices/{Tissue}.{ri_group}.mat.gz', Tissue=tissue_sub_list, ri_group = ['snmf_and_gencode.retained_introns',
               'snmf_only.retained_introns', 'gencode_only.retained_introns', 'appris_introns', 'snmf.flanking_exons_1', 'gencode.flanking_exons_1',
               'appris_flanking_exons_1', 'snmf.flanking_exons_2', 'gencode.flanking_exons_2', 'appris_flanking_exons_2', 'snmf.retained_introns'])







rule GetInputForIntronMetaplot_split:
    input:
        snmf = 'CoveragePlots/bed_files/introns/snmf.retained_introns.bed.gz',
        gencode = 'CoveragePlots/bed_files/introns/gencode.retained_introns.bed.gz',
        gencode_exons = 'Annotations/gencode.v44.primary_assembly.exons.bed.gz', #'Annotations/gencode.v44.primary_assembly.exons.sorted.bed.gz',
        snmf_exons = 'ebpmf_models/filtered/snmf_10/tables/snmf.merged_isoforms.exons.sorted.bed.gz',
        appris = 'CoveragePlots/bed_files/introns/gencode.appris_introns.bed.gz',
        vastdb = 'CoveragePlots/bed_files/introns/vastdb.retained_introns.bed.gz'
    output:
        snmf_only_temp = temp('CoveragePlots/bed_files/introns_split/snmf_only_temp.retained_introns.bed'),
        snmf_only = 'CoveragePlots/bed_files/introns_split/snmf_only.retained_introns.bed',
        gencode_only_temp = temp('CoveragePlots/bed_files/introns_split/gencode_only_temp.retained_introns.bed'),
        gencode_only = 'CoveragePlots/bed_files/introns_split/gencode_only.retained_introns.bed',
        snmf_gencode = 'CoveragePlots/bed_files/introns_split/snmf_and_gencode.retained_introns.bed',
        appris = 'CoveragePlots/bed_files/introns_split/appris_introns.bed',
        vastdb_only_temp = temp('CoveragePlots/bed_files/introns_split/vastdb_only_temp.retained_introns.bed'),
        vastdb_only = 'CoveragePlots/bed_files/introns_split/vastdb_only.retained_introns.bed',
        snmf_vastdb = 'CoveragePlots/bed_files/introns_split/snmf_and_vastdb.retained_introns.bed',
        gencode_vastdb = 'CoveragePlots/bed_files/introns_split/gencode_and_vastdb.retained_introns.bed',
        snmf_gencode_vastdb = 'CoveragePlots/bed_files/introns_split/snmf_and_gencode_and_vastdb.retained_introns.bed'
    log:
        "/scratch/midway3/cnajar/logs/bigwigs/get_introns_bed_for_metaplots_split.log"
    resources:
        mem_mb = 12000
    shell:
        """
        (bedtools intersect -v -a {input.snmf} -b {input.gencode} | bedtools intersect -v -a - -b {input.vastdb} | bedtools sort -i - > {output.snmf_only_temp}) &>> {log};
        (bedtools intersect -v -f 0.01 -a {output.snmf_only_temp} -b {input.gencode_exons} > {output.snmf_only}) &>> {log};
        (bedtools intersect -v -a {input.gencode} -b {input.snmf} | bedtools intersect -v -a - -b {input.vastdb}  | bedtools sort -i - > {output.gencode_only_temp}) &>> {log};
        (bedtools intersect -v -f 0.01 -a {output.gencode_only_temp} -b {input.snmf_exons} > {output.gencode_only}) &>> {log};
        (bedtools intersect -v -a {input.vastdb} -b {input.snmf} | bedtools intersect -v -a - -b {input.gencode}  | bedtools sort -i - > {output.vastdb_only_temp}) &>> {log};
        (bedtools intersect -v -f 0.01 -a {output.vastdb_only_temp} -b {input.snmf_exons} | bedtools intersect -v -f 0.01 -a - -b {input.gencode_exons} > {output.vastdb_only}) &>> {log};
        (bedtools intersect -s -F 1 -f 1 -u -a {input.snmf} -b {input.gencode} | bedtools intersect -v -a - -b {input.vastdb} | bedtools sort -i - > {output.snmf_gencode}) &>> {log};
        (bedtools intersect -s -F 1 -f 1 -u -a {input.snmf} -b {input.vastdb} | bedtools intersect -v -f 0.01 -a - -b {input.gencode_exons} | bedtools sort -i - > {output.snmf_vastdb}) &>> {log};
        (bedtools intersect -s -F 1 -f 1 -u -a {input.gencode} -b {input.vastdb} | bedtools intersect -v -f 0.01 -a - -b {input.snmf_exons} | bedtools sort -i - > {output.gencode_vastdb}) &>> {log};
        (bedtools intersect -s -F 1 -f 1 -u -a {input.snmf} -b {input.gencode} | bedtools intersect -s -F 1 -f 1 -u -a - -b {input.vastdb} | bedtools sort -i - > {output.snmf_gencode_vastdb}) &>> {log};
        (bedtools intersect -v -f 0.01 -F 0.01 -a {input.appris} -b {input.snmf_exons} | bedtools intersect -v -f 0.01 -F 0.01 -a - -b {input.gencode_exons} | bedtools intersect -v -f 0.01 -F 0.01 -a - -b {input.vastdb} | bedtools sort -i - > {output.appris}) &>> {log}
        """

use rule ComputeMatrixForCoverageMetaplots_sep as ComputeMatrixForCoverageMetaplots_sep_split with:
    input:
        bigwigs = "CoveragePlots/bigwigs/{Tissue}.bw",
        bed = "CoveragePlots/bed_files/introns_split/{ri_group}.bed"
    output:
        mat = "CoveragePlots/matrices_split/{Tissue}.{ri_group}.mat.gz"
    wildcard_constraints:
        Tissue = "|".join(tissue_sub_list),
        ri_group = '|'.join(['snmf_and_gencode.retained_introns', 'snmf_only.retained_introns', 'gencode_only.retained_introns', 'appris_introns',
'vastdb_only.retained_introns', 'snmf_and_vastdb.retained_introns', 'gencode_and_vastdb.retained_introns', 'snmf_and_gencode_and_vastdb.retained_introns'])
    log:
        "logs/Metaplots/ComputeMatrix_sep_split.{Tissue}.{ri_group}.log"

#rule ComputeMatrixForCoverageMetaplots_sep_split:
#    input:
#        bigwigs = "CoveragePlots/bigwigs/{Tissue}.bw",
#        bed = "CoveragePlots/bed_files/introns_split/{ri_group}.bed" 
#    output:
#        mat = "CoveragePlots/matrices_split/{Tissue}.{ri_group}.mat.gz"
#    wildcard_constraints:
#        Tissue = "|".join(tissue_sub_list),
#        ri_group = '|'.join(['snmf_and_gencode.retained_introns', 'snmf_only.retained_introns', 'gencode_only.retained_introns', 'appris_introns',
#'vastdb_only.retained_introns', 'snmf_and_vastdb.retained_introns', 'gencode_and_vastdb.retained_introns', 'snmf_and_gencode_and_vastdb.retained_introns'])
#    log:
#        "logs/Metaplots/ComputeMatrix.{Tissue}.{ri_group}.log"
#    resources:
#        mem_mb = 12000
#    shell:
#        """
#        computeMatrix scale-regions -S {input.bigwigs} -R {input.bed} -b 150 -a 150 --startLabel "5' splice site" --endLabel "3' splice site" -o {output} &> #{log}
#        """


rule collect_metaplot_matrices_split:
    input:
        expand('CoveragePlots/matrices_split/{Tissue}.{ri_group}.mat.gz', Tissue=tissue_sub_list, ri_group = ['snmf_and_gencode.retained_introns', 'snmf_only.retained_introns', 'gencode_only.retained_introns', 'appris_introns', 'vastdb_only.retained_introns', 'snmf_and_vastdb.retained_introns', 'gencode_and_vastdb.retained_introns', 'snmf_and_gencode_and_vastdb.retained_introns'])




rule GetCassetteExonsBeds:
    input:
        'Annotations/gencode.v44.primary_assembly.exons.sorted.bed.gz',
        'ebpmf_models/filtered/snmf_10/tables/snmf.merged_isoforms.exons.sorted.bed.gz',
        'ebpmf_models/filtered/snmf_10/tables/annotated.snmf.merged_isoforms.tab.gz',
        'coverage/GTEx_overlapping_regions.bed.gz'
    output:
        expand('CoveragePlots/bed_files/cassette_exons/{annotation}/{capture}.{region}.bed', 
                annotation = ['gencode', 'vastdb'], capture = ['exon_only', 'cassette_exon', 'intron_only'], 
                region = ['exon', 'flanking1', 'flanking2']),
        expand('CoveragePlots/bed_files/cassette_exons/snmf/{annotation_match}.{region}.bed', 
                annotation_match = ['gencode', 'vastdb', 'junctions', 'unmatched'], region = ['exon', 'flanking1', 'flanking2']),
        'ebpmf_models/filtered/snmf_10/tables/cassette_exons.bed.gz'
    log:
        "/scratch/midway3/cnajar/logs/bigwigs/get_cassette_exons_bed.log"
    resources:
        mem_mb = 52000
    shell:
        """
        python scripts/get_cassette_exons.py &> {log}
        """

rule CollectCassetteExons:
    input:
        expand('CoveragePlots/bed_files/cassette_exons/{annotation}/{capture}.{region}.bed', 
                annotation = ['gencode', 'vastdb'], capture = ['exon_only', 'cassette_exon', 'intron_only'], 
                region = ['exon', 'flanking1', 'flanking2']),
        expand('CoveragePlots/bed_files/cassette_exons/snmf/{annotation_match}.{region}.bed', 
                annotation_match = ['gencode', 'vastdb', 'junctions', 'unmatched'], region = ['exon', 'flanking1', 'flanking2']),
        'ebpmf_models/filtered/snmf_10/tables/cassette_exons.bed.gz'

rule ComputeMatrixForCoverageMetaplots_CassetteExons:
    input:
        bigwigs = "CoveragePlots/bigwigs/{Tissue}.bw",
        bed = 'CoveragePlots/bed_files/cassette_exons/{annotation}/{capture}.{region}.bed',
    output:
        mat = "CoveragePlots/matrices/cassette_exons/{Tissue}.{annotation}.{capture}.{region}.mat.gz"
    wildcard_constraints:
        Tissue = "|".join(tissue_sub_list),
        annotation = 'gencode|vastdb',
        capture = 'exon_only|cassette_exon|intron_only',
        region = 'exon|flanking1|flanking2'
    log:
        "logs/Metaplots/ComputeMatrix.CE.{Tissue}.{annotation}.{capture}.{region}.log"
    resources:
        mem_mb = 12000
    params:
        '--regionBodyLength 300 --binSize 5'
    shell:
        """
        computeMatrix scale-regions -S {input.bigwigs} -R {input.bed} {params} -b 150 -a 150 --startLabel "5' splice site" --endLabel "3' splice site" -o {output} &> {log}
        """

use rule ComputeMatrixForCoverageMetaplots_CassetteExons as ComputeMatrixForCoverageMetaplots_CassetteExons_snmf with:
    input:
        bigwigs = "CoveragePlots/bigwigs/{Tissue}.bw",
        bed = 'CoveragePlots/bed_files/cassette_exons/snmf/{annotation_match}.{region}.bed'
    output:
        mat = "CoveragePlots/matrices/cassette_exons/{Tissue}.snmf.{annotation_match}.{region}.mat.gz"
    wildcard_constraints:
        Tissue = "|".join(tissue_sub_list),
        annotation_match = 'gencode|vastdb|junctions|unmatched',
        region = 'exon|flanking1|flanking2'
    log:
        "logs/Metaplots/ComputeMatrix.CE.{Tissue}.snmf.{annotation_match}.{region}.log"

rule CollectCassetteExons_metaplots:
    input:
        expand("CoveragePlots/matrices/cassette_exons/{Tissue}.{annotation}.{capture}.{region}.mat.gz", Tissue = tissue_sub_list, 
                annotation = ['gencode', 'vastdb'], capture = ['exon_only', 'cassette_exon', 'intron_only'], 
                region = ['exon', 'flanking1', 'flanking2']),
        expand("CoveragePlots/matrices/cassette_exons/{Tissue}.snmf.{annotation_match}.{region}.mat.gz", Tissue = tissue_sub_list,
                annotation_match = ['gencode', 'vastdb', 'junctions', 'unmatched'], region = ['exon', 'flanking1', 'flanking2']),

rule GetCoverageGaps:
    output:
        'coverage/GTEx_overlapping_regions.bed.gz'
    log:
        'logs/coverage_gaps.log'
    resources:
        mem_mb = 8000
    shell:
        """
        python scripts/get_gaps_in_gtex_models.py &> {log}
        """
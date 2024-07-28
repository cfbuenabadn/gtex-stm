qqnorm_output = []
for tissue_ in tissue_sub_list:
    for K in [3, 5, 10]:
        qqnorm_output.append(f'QTLs/GTEx_10/{tissue_}/{{annotation}}_{str(K)}.qqnorm.bed.gz')

def GetQQNormInput(wildcards):
    if wildcards.annotation == 'snmf':
        return expand('ebpmf_models/filtered/snmf_{K}/tables/EL.bed.gz', K=[3, 5, 10])
    elif wildcards.annotation == 'transcripts':
        return expand('ebpmf_models/filtered/snmf_{K}/tables/transcript.merged_isoforms.EL.bed.gz', K=[3, 5, 10])

rule QQNormTissueTables_GTEx_10:
    input:
        GetQQNormInput #expand('ebpmf_models/filtered/snmf_{K}/tables/EL.bed.gz', K=[3, 5, 10]),
    output:
        qqnorm_output
    log:
        "/scratch/midway3/cnajar/logs/QTLs/GTEx_10/QQNormsNMF_Tables.{annotation}.log"
    resources:
        mem_mb = 24000
    wildcard_constraints:
        annotation = 'snmf|transcripts'
    shell:
        """
        python scripts/QQnorm_sNMF_tables.py {wildcards.annotation} &> {log}
        """

rule SortQTLtoolsPhenotypeTable:
    input:
        'QTLs/GTEx_10/{tissue_id}/{annotation}_{K}.qqnorm.bed.gz'
    output:
        bed = 'QTLs/GTEx_10/{tissue_id}/{annotation}_{K}.sorted.qqnorm.bed.gz',
        tbi = 'QTLs/GTEx_10/{tissue_id}/{annotation}_{K}.sorted.qqnorm.bed.gz.tbi'
    log:
        "/scratch/midway3/cnajar/logs/QTLs/SortQTLtoolsPhenotypeTable/{annotation}.{tissue_id}.{K}.log"
    resources:
        mem_mb = 24000
    wildcard_constraints:
        tissue_id = '|'.join(tissue_list),
        K = '2|3|4|5|10',
        annotation = 'snmf|transcripts'
    shell:
        """
        (bedtools sort -header -i {input} | bgzip /dev/stdin -c > {output.bed}) &> {log}
        (tabix -p bed {output.bed}) &>> {log}
        """


rule PhenotypePCs:
    """
    QTLtools format expression PCs as covariates
    including the number of PCs that explain more
    variance then when the phenotype table is
    permuted
    """
    input:
        'QTLs/GTEx_10/{tissue_id}/{annotation}_{K}.sorted.qqnorm.bed.gz',
    output:
        'QTLs/GTEx_10/{tissue_id}/{annotation}_{K}.sorted.qqnorm.bed.pca',
    log:
        "/scratch/midway3/cnajar/logs/QTLs/PhenotypePCs/{tissue_id}.{annotation}_{K}.log"
    resources:
        mem_mb = 24000
    wildcard_constraints:
        annotation = 'snmf|transcripts',
        tissue_id = '|'.join(tissue_list),
        K = '2|3|4|5|10'
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/PermuteAndPCA.R {input} {output} &> {log}
        """

rule collect_PCA:
    input:
        expand('QTLs/GTEx_10/{tissue_id}/{annotation}_{K}.sorted.qqnorm.bed.pca', annotation = ['snmf', 'transcripts'], 
               tissue_id = tissue_sub_list, K = [2, 3, 4, 5, 10]),
       

N_PermutationChunks = 100
ChunkNumbers = range(0, 1+N_PermutationChunks) 


def GetQTLtoolsPassFlags(wildcards):
    if wildcards.Pass in ["GroupedPermutationPass", "PermutationPass"]:
        return "--permute 1000"
    elif wildcards.Pass in ["NominalPass", "GroupedNominalPass"]:
        return "--nominal 1"

# In case I need to exclude outliers
def GetExcludeFile(wildcards):
    if wildcards.SexGroup == '.female':
        return "--exclude-samples config/female_participants.txt"
    elif wildcards.SexGroup == '.male':
        return "--exclude-samples config/male_participants.txt"
    else:
        return ""

def GetQTLtoolsWindowFlag(wildcards):
    if wildcards.FeatureCoordinatesRedefinedFor == '.ForGWASColoc':
        return "--window 1"
    else:
        "--window 10000"


rule MakePhenotypeTableToColocFeaturesWithGWASLoci:
    """
    To colocalize molQTLs with gwas we need summary stats (betas and se) for
    every snp in the same window as a gwas locus. Here, I will take a phenotype
    table, find the features that intersect the gwas locus (1MB window centered
    on lead snp, defined in the input file), and output a phenotype table with
    every intersection with coordinates redefined to be the 1MB gwas locus
    window. Running QTLtools on this will help me get the necessary summary
    stats in those windows
    """
    input:
        bed = "QTLs/GTEx_10/{tissue}/{annotation}_{K}.sorted.qqnorm.bed.gz",
        fai = "Annotations/GRCh38.primary_assembly.genome.fa.fai",
        loci = "gwas/LeadSnpWindows.bed"
    output:
        bed = "QTLs/GTEx_10/{tissue}/{annotation}_{K}.ForGWASColoc.sorted.qqnorm.bed.gz",
        tbi = "QTLs/GTEx_10/{tissue}/{annotation}_{K}.ForGWASColoc.sorted.qqnorm.bed.gz.tbi"
    resources:
        mem_mb = 48000
    wildcard_constraints:
        annotation = 'snmf|transcripts',
    log:
        "/scratch/midway3/cnajar/logs/MakePhenotypeTableToColocFeaturesWithGWASLoci/{annotation}.{tissue}.{K}.log"
    shell:
        """
        (cat <(zcat {input.bed} | head -1) <(  bedtools intersect  -wo -a {input.bed} -b {input.loci} -sorted | awk -F'\\t' -v OFS='\\t' '{{$4=$4":"$(NF-1); $5=$(NF-1); $2=$(NF-3); $3=$(NF-2); print $0}}' | rev | cut -f 6- | rev ) | tr ' ' '\\t' | bedtools sort -i - -header | bgzip /dev/stdin -c > {output.bed} ) &> {log}
        tabix -p bed {output.bed}
        """

def GetOtherFlags(wildcards):
    if wildcards.Pass in ["GroupedPermutationPass", "GroupedNominalPass"]:
        return "--grp-best"
    else:
        return ""

rule SelectSignificantTraitsForGWASColoc:
    input:
        bed = 'QTLs/GTEx_10/{tissue_id}/{annotation}_{K}{FeatureCoordinatesRedefinedFor}.sorted.qqnorm.bed.gz',
        qtls = 'QTLs/GTEx_10/{tissue_id}/{annotation}_{K}{FeatureCoordinatesRedefinedFor}.PermutationPass.FDR_Added.txt.gz'
    output:
        bed = temp('QTLs/GTEx_10/{tissue_id}/{annotation}_{K}{FeatureCoordinatesRedefinedFor}.sorted.significant_only.qqnorm.bed'),
        bgz = 'QTLs/GTEx_10/{tissue_id}/{annotation}_{K}{FeatureCoordinatesRedefinedFor}.sorted.significant_only.qqnorm.bed.gz',
        tbi = 'QTLs/GTEx_10/{tissue_id}/{annotation}_{K}{FeatureCoordinatesRedefinedFor}.sorted.significant_only.qqnorm.bed.gz.tbi'
    log:
        '/scratch/midway3/cnajar/logs/QTLs/get_significant_traits.ForGWASColoc.{tissue_id}.{annotation}_{K}{FeatureCoordinatesRedefinedFor}.log'
    wildcard_constraints:
        annotation = 'snmf|transcripts',
        K = '2|3|4|5|10',
        FeatureCoordinatesRedefinedFor = '.ForGWASColoc'
    shell:
        """
        python scripts/filter_permutation_pass_forGWAScoloc.py {input.bed} {input.qtls} {output.bed} &> {log};
        (bgzip {output.bed} -c > {output.bgz}) &>> {log}
        tabix -p bed {output.bgz} &>> {log}
        """

def GetQTLToolsBedInput(wildcards):
    if (wildcards.Pass == 'NominalPass') and (wildcards.FeatureCoordinatesRedefinedFor == '.ForGWASColoc'):
        return 'QTLs/GTEx_10/{tissue_id}/{annotation}_{K}{FeatureCoordinatesRedefinedFor}.sorted.significant_only.qqnorm.bed.gz'
    else:
        return 'QTLs/GTEx_10/{tissue_id}/{annotation}_{K}{FeatureCoordinatesRedefinedFor}.sorted.qqnorm.bed.gz'

def GetQTLToolsTbiInput(wildcards):
    if (wildcards.Pass == 'NominalPass') and (wildcards.FeatureCoordinatesRedefinedFor == '.ForGWASColoc'):
        return 'QTLs/GTEx_10/{tissue_id}/{annotation}_{K}{FeatureCoordinatesRedefinedFor}.sorted.significant_only.qqnorm.bed.gz.tbi'
    else:
        return 'QTLs/GTEx_10/{tissue_id}/{annotation}_{K}{FeatureCoordinatesRedefinedFor}.sorted.qqnorm.bed.gz.tbi'

rule QTLtools_generalized:
    input:
        vcf = "/project2/yangili1/cfbuenabadn/ChromatinSplicingQTLs/code/GTEx/data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.vcf.gz",
        tbi = "/project2/yangili1/cfbuenabadn/ChromatinSplicingQTLs/code/GTEx/data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.vcf.gz.tbi",
        bed = GetQTLToolsBedInput, #'QTLs/GTEx_10/{tissue_id}/{annotation}_{K}{FeatureCoordinatesRedefinedFor}.sorted.qqnorm.bed.gz',
        bed_tbi = GetQTLToolsTbiInput, #'QTLs/GTEx_10/{tissue_id}/{annotation}_{K}{FeatureCoordinatesRedefinedFor}.sorted.qqnorm.bed.gz.tbi',
        cov = 'QTLs/GTEx_10/{tissue_id}/{annotation}_{K}.sorted.qqnorm.bed.pca'
    output:
        "QTLs/GTEx_10/{tissue_id}/{annotation}_{K}{SexGroup}{FeatureCoordinatesRedefinedFor}.{Pass}.Chunks/{QTLTools_chunk_n}.txt"
    log:
        "/scratch/midway3/cnajar/logs/QTLs/{tissue_id}/{annotation}_{K}{SexGroup}{FeatureCoordinatesRedefinedFor}.{Pass}.Chunks/{QTLTools_chunk_n}.log"
    resources:
        mem_mb = much_more_mem_after_first_attempt
    envmodules:
        "gsl/2.5"
    params:
        WindowFlag = GetQTLtoolsWindowFlag,#"--window 10000",
        OtherFlags = GetOtherFlags,#"--grp-best",
        PassFlags = GetQTLtoolsPassFlags,
        ExcFlag = GetExcludeFile,
        NChunks = N_PermutationChunks
    wildcard_constraints:
        annotation = 'snmf|transcripts',
        n = "|".join(str(i) for i in ChunkNumbers),
        Pass = "GroupedPermutationPass|PermutationPass|GroupedNominalPass|NominalPass",
        SexGroup = "|.female|.male",
        K = '2|3|4|5|10',
        FeatureCoordinatesRedefinedFor = '|.ForGWASColoc'
    shell:
        """
        {config[QTLtools]} cis --std-err --chunk {wildcards.QTLTools_chunk_n} {params.NChunks} --vcf {input.vcf} --bed {input.bed} --cov {input.cov} --out {output} {params.OtherFlags} {params.WindowFlag} {params.PassFlags} {params.ExcFlag} &> {log}
        if [ ! -f {output} ]
        then
            touch {output}
        fi
        """



rule Gather_QTLtools_cis_pass:
    input:
        expand( "QTLs/GTEx_10/{{tissue_id}}/{{annotation}}_{{K}}{{SexGroup}}{{FeatureCoordinatesRedefinedFor}}.{{Pass}}.Chunks/{QTLTools_chunk_n}.txt", QTLTools_chunk_n=ChunkNumbers )
    output:
        "QTLs/GTEx_10/{tissue_id}/{annotation}_{K}{SexGroup}{FeatureCoordinatesRedefinedFor}.{Pass}.txt.gz"
    log:
        "/scratch/midway3/cnajar/logs/Gather_QTLtools_cis_pass/{tissue_id}{SexGroup}{FeatureCoordinatesRedefinedFor}.{annotation}_{K}.{Pass}.log"
    wildcard_constraints:
        annotation = 'snmf|transcripts',
        Pass = "GroupedPermutationPass|PermutationPass|GroupedNominalPass|NominalPass",
        SexGroup = "|.female|.male",
        FeatureCoordinatesRedefinedFor = '|.ForGWASColoc'
    resources:
        mem_mb = much_more_mem_after_first_attempt
    shell:
        """
        (cat {input} | gzip - > {output}) &> {log}
        """
        

rule AddQValueToPermutationPass:
    input:
        "QTLs/GTEx_10/{tissue_id}/{annotation}_{K}{SexGroup}{FeatureCoordinatesRedefinedFor}.{Pass}.txt.gz"
    output:
        table = "QTLs/GTEx_10/{tissue_id}/{annotation}_{K}{SexGroup}{FeatureCoordinatesRedefinedFor}.{Pass}.FDR_Added.txt.gz",
    log:
        "/scratch/midway3/cnajar/logs/AddQValueToPermutationPass/{tissue_id}/{SexGroup}{FeatureCoordinatesRedefinedFor}.{annotation}_{K}.{Pass}.txt.gz.log"
    conda:
        "../envs/r_essentials.yml"
    wildcard_constraints:
        annotation = 'snmf|transcripts',
        SexGroup = "|.female|.male",
        FeatureCoordinatesRedefinedFor = '|.ForGWASColoc',
        Pass = "PermutationPass|GroupedPermutationPass"
    priority:
        10
    shell:
        """
        Rscript scripts/AddQvalueToQTLtoolsOutput.R {input} {output.table} {wildcards.Pass} &> {log}
        """


rule CollectQTLs:
    input:
        expand("QTLs/GTEx_10/{tissue_id}/{annotation}_{K}{SexGroup}.{Pass}.txt.gz",
               annotation = ['snmf'],
               Pass = ['GroupedPermutationPass.FDR_Added', 'GroupedNominalPass'],
               tissue_id = tissue_sub_list,
               SexGroup = [''],# '.female', '.male'],
               K = [3, 10]),
        expand("QTLs/GTEx_10/{tissue_id}/{annotation}_{K}{SexGroup}.GroupedNominalPass.txt.tabix.gz",
               annotation = ['snmf'],
               tissue_id = tissue_sub_list,
               SexGroup = ['', '.female', '.male'],
               K = [3, 10]),
        expand("QTLs/GTEx_10/{tissue_id}/{annotation}_{K}{SexGroup}.{Pass}.txt.gz",
               annotation = ['transcripts'],
               Pass = ['GroupedPermutationPass.FDR_Added', 'GroupedNominalPass'],
               tissue_id = tissue_sub_list,
               SexGroup = [''],# '.female', '.male'],
               K = [10]),
        expand("QTLs/GTEx_10/{tissue_id}/{annotation}_{K}{SexGroup}.GroupedNominalPass.txt.tabix.gz",
               annotation = ['transcripts'],
               tissue_id = tissue_sub_list,
               SexGroup = ['', '.female', '.male'],
               K = [10])

rule tabixNominalPassQTLResults:
    """
    Convert QTLtools output to tab delimited bgzipped and tabix indexed files
    for easy access with tabix
    """
    input:
        "QTLs/GTEx_10/{tissue_id}/{annotation}_{K}{SexGroup}.GroupedNominalPass.txt.gz"
    wildcard_constraints:
        annotation = 'snmf|transcripts',
        tissue_id = "|".join(tissue_sub_list),
        SexGroup = "|.male|.female",
        K = '2|3|4|5|10',
        #FeatureCoordinatesRedefinedFor = '',
        #Pass = 'GroupedNominalPass'
    params:
        sort_temp = '-T ' + config['scratch'][:-1]
        # sort_temp = ""
    output:
        txt = "QTLs/GTEx_10/{tissue_id}/{annotation}_{K}{SexGroup}.GroupedNominalPass.txt.tabix.gz",
        tbi = "QTLs/GTEx_10/{tissue_id}/{annotation}_{K}{SexGroup}.GroupedNominalPass.txt.tabix.gz.tbi"
    resources:
        mem_mb = much_more_mem_after_first_attempt
    log:
        "/scratch/midway3/cnajar/logs/tabixNominalPassQTLResults/{annotation}.{tissue_id}{SexGroup}.snmf_{K}.GroupedNominalPass.log"
    shadow: "shallow"
    shell:
        """
        (cat <(zcat {input} | head -1 | perl -p -e 'printf("#") if $. ==1; s/ /\\t/g') <(zcat {input} | awk 'NR>1' |  perl -p -e 's/ /\\t/g' | sort {params.sort_temp} -k11,11 -k12,12n  ) | bgzip /dev/stdin -c > {output.txt}) &> {log}
        tabix -b 12 -e12 -s11 {output.txt} &>> {log}
        """

rule tabixNominalPassQTLResultsForGWASColoc:
    """
    Convert QTLtools output to tab delimited bgzipped and tabix indexed files
    for easy access with tabix
    """
    input:
        "QTLs/GTEx_10/{tissue_id}/{annotation}_{K}{SexGroup}{FeatureCoordinatesRedefinedFor}.{Pass}.txt.gz"
    wildcard_constraints:
        annotation = 'snmf|transcripts',
        tissue_id = "|".join(tissue_sub_list),
        SexGroup = "|.male|.female",
        K = '2|3|4|5|10',
        FeatureCoordinatesRedefinedFor = '.ForGWASColoc',
        Pass = 'NominalPass'
    params:
        sort_temp = '-T ' + config['scratch'][:-1]
        # sort_temp = ""
    output:
        txt = "QTLs/GTEx_10/{tissue_id}/{annotation}_{K}{SexGroup}{FeatureCoordinatesRedefinedFor}.{Pass}.txt.tabix.gz",
        tbi = "QTLs/GTEx_10/{tissue_id}/{annotation}_{K}{SexGroup}{FeatureCoordinatesRedefinedFor}.{Pass}.txt.tabix.gz.tbi"
    resources:
        mem_mb = much_more_mem_after_first_attempt
    log:
        "/scratch/midway3/cnajar/logs/tabixNominalPassQTLResults/{tissue_id}{SexGroup}{FeatureCoordinatesRedefinedFor}.{annotation}_{K}.{Pass}.log"
    shadow: "shallow"
    shell:
        """
        (cat <(zcat {input} | head -1 | perl -p -e 'printf("#") if $. ==1; s/ /\\t/g') <(zcat {input} | awk 'NR>1' |  perl -p -e 's/ /\\t/g' | sort {params.sort_temp} -k9,9 -k10,10n  ) | bgzip /dev/stdin -c > {output.txt}) &> {log}
        tabix -b 10 -e10 -s9 {output.txt} &>> {log}
        """



rule CollectQTLsForGWASColoc:
    input:
        expand("QTLs/GTEx_10/{tissue_id}/{annotation}_{K}.ForGWASColoc.PermutationPass.FDR_Added.txt.gz", 
               annotation = ['snmf'],
               tissue_id = tissue_sub_list,
               K = [3, 10]),
        expand("QTLs/GTEx_10/{tissue_id}/{annotation}_{K}.ForGWASColoc.NominalPass.txt.tabix.gz",
               annotation = ['snmf'],
               tissue_id = tissue_sub_list,
               K = [3, 10])



tissue_sub_list_no_heart = [x for x in tissue_sub_list if "Heart" not in x]

rule CollectTranscriptQTLs:
    input:
        expand("QTLs/GTEx_10/{tissue_id}/transcripts_10.{Pass}.txt.gz",
               Pass = ['GroupedPermutationPass.FDR_Added', 'ForGWASColoc.PermutationPass.FDR_Added'],
               tissue_id = tissue_sub_list),
         expand("QTLs/GTEx_10/{tissue_id}/transcripts_10.ForGWASColoc.NominalPass.txt.tabix.gz",
               tissue_id = tissue_sub_list),
         #expand("QTLs/GTEx_10/{tissue_id}/snmf_3.{Pass}.txt.gz",
         #      Pass = ['GroupedPermutationPass.FDR_Added', 'ForGWASColoc.PermutationPass.FDR_Added'],
         #      tissue_id = tissue_sub_list),
         #expand("QTLs/GTEx_10/{tissue_id}/snmf_3.ForGWASColoc.NominalPass.txt.tabix.gz",
         #      tissue_id = tissue_sub_list)









rule sQTLsMakeInput:
    input:
        expand('GTEx_Analysis_v8_sQTL_phenotype_matrices/{tissue_id}.v8.leafcutter_phenotypes.bed.gz', tissue_id = tissue_sub_list)
    output:
        expand('QTLs/sQTLs/{tissue_id}/leafcutter.qqnorm.bed.gz', tissue_id = tissue_sub_list)
    log:
        '/scratch/midway3/cnajar/logs/make.sqtl.input.log'
    shell:
        """
        python scripts/MakeQTLToolsSplicingInputs.py &> {log}
        """

use rule SortQTLtoolsPhenotypeTable as SortQTLtoolsPhenotypeTable_sQTL with:
    input:
        'QTLs/sQTLs/{tissue_id}/leafcutter.qqnorm.bed.gz'
    output:
        bed = 'QTLs/sQTLs/{tissue_id}/leafcutter.sorted.qqnorm.bed.gz',
        tbi = 'QTLs/sQTLs/{tissue_id}/leafcutter.sorted.qqnorm.bed.gz.tbi'
    log:
        "/scratch/midway3/cnajar/logs/QTLs/SortQTLtoolsPhenotypeTable/sqtl.{tissue_id}.log"
    resources:
        mem_mb = 32000

use rule PhenotypePCs as PhenotypePCs_sQTLs with:
    input:
        'QTLs/sQTLs/{tissue_id}/leafcutter.sorted.qqnorm.bed.gz',
    output:
        'QTLs/sQTLs/{tissue_id}/leafcutter.sorted.qqnorm.bed.pca',
    log:
        "/scratch/midway3/cnajar/logs/QTLs/PhenotypePCs/sqtls.{tissue_id}.log"

use rule QTLtools_generalized as sQTLs with:
    input:
        vcf = "/project2/yangili1/cfbuenabadn/ChromatinSplicingQTLs/code/GTEx/data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.vcf.gz",
        tbi = "/project2/yangili1/cfbuenabadn/ChromatinSplicingQTLs/code/GTEx/data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.vcf.gz.tbi",
        bed = 'QTLs/sQTLs/{tissue_id}/leafcutter.sorted.qqnorm.bed.gz',
        bed_tbi = 'QTLs/sQTLs/{tissue_id}/leafcutter.sorted.qqnorm.bed.gz.tbi',
        cov = 'QTLs/sQTLs/{tissue_id}/leafcutter.sorted.qqnorm.bed.pca'
    output:
        temp("QTLs/sQTLs/{tissue_id}/leafcutter.{Pass}.Chunks/{QTLTools_chunk_n}.txt")
    log:
        "/scratch/midway3/cnajar/logs/QTLs/sQTLs/{tissue_id}/leafcutter.{Pass}.Chunks/{QTLTools_chunk_n}.log"
    params:
        WindowFlag = "--window 10000",
        OtherFlags = "--grp-best",
        PassFlags = GetQTLtoolsPassFlags,
        ExcFlag = "",
        NChunks = N_PermutationChunks
    wildcard_constraints:
        n = "|".join(str(i) for i in ChunkNumbers),
        Pass = "GroupedPermutationPass|GroupedNominalPass"

use rule Gather_QTLtools_cis_pass as Gather_QTLtools_cis_pass_sqtls with:
    input:
        expand( "QTLs/sQTLs/{{tissue_id}}/leafcutter.{{Pass}}.Chunks/{QTLTools_chunk_n}.txt", QTLTools_chunk_n=ChunkNumbers )
    output:
        "QTLs/sQTLs/{tissue_id}/leafcutter.{Pass}.txt.gz"
    log:
        "/scratch/midway3/cnajar/logs/Gather_QTLtools_cis_pass/leafcutter.{tissue_id}.{Pass}.log"
    wildcard_constraints:

use rule AddQValueToPermutationPass as AddQValueToPermutationPass_sqtls with:
    input:
        "QTLs/sQTLs/{tissue_id}/leafcutter.{Pass}.txt.gz"
    output:
        table = "QTLs/sQTLs/{tissue_id}/leafcutter.{Pass}.FDR_Added.txt.gz",
    log:
        "/scratch/midway3/cnajar/logs/AddQValueToPermutationPass/{tissue_id}/leafcutter.{Pass}.txt.gz.log"

use rule tabixNominalPassQTLResults as tabixNominalPassQTLResults_sqtls with:
    input:
        "QTLs/sQTLs/{tissue_id}/leafcutter.GroupedNominalPass.txt.gz"
    params:
        sort_temp = '-T ' + config['scratch'][:-1]
    output:
        txt = "QTLs/sQTLs/{tissue_id}/leafcutter.GroupedNominalPass.txt.tabix.gz",
        tbi = "QTLs/sQTLs/{tissue_id}/leafcutter.GroupedNominalPass.txt.tabix.gz.tbi"
    log:
        "/scratch/midway3/cnajar/logs/tabixNominalPassQTLResults/leafcutter.{tissue_id}.GroupedNominalPass.log"


rule Collect_sQTLs:
    input:
        expand("QTLs/sQTLs/{tissue_id}/leafcutter.{Pass}.txt.gz", 
               Pass = ['GroupedPermutationPass.FDR_Added', 'GroupedNominalPass'],
               tissue_id = tissue_sub_list),
        expand("QTLs/sQTLs/{tissue_id}/leafcutter.GroupedNominalPass.txt.tabix.gz",
               tissue_id = tissue_sub_list)
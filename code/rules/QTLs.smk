

rule QQNormTissueTables:
    input:
        expand('ebpmf_models/filtered/snmf_{k}/tables/EL.bed.gz', k=[2, 3, 4, 5, 10]),
        'ebpmf_models/single_tissue/tables/BA9.EL.bed.gz',
        'ebpmf_models/single_tissue/tables/MS.EL.bed.gz',
        'ebpmf_models/single_tissue/tables/WB.EL.bed.gz',
        'ebpmf_models/single_tissue/tables/Skin.EL.bed.gz'
    output:
        qqnorm_output
    log:
        "logs/QTLs/QQNormsNMF_Tables.log"
    resources:
        mem_mb = 24000
    shell:
        """
        python scripts/QQnorm_sNMF_tables.py &> {log}
        """

rule SortQTLtoolsPhenotypeTable:
    input:
        'QTLs/{tissue_id}/{MultiOrSingle}.snmf_{K}.qqnorm.bed.gz'
    output:
        bed = 'QTLs/{tissue_id}/{MultiOrSingle}.snmf_{K}.sorted.qqnorm.bed.gz',
        tbi = 'QTLs/{tissue_id}/{MultiOrSingle}.snmf_{K}.sorted.qqnorm.bed.gz.tbi'
    log:
        "logs/QTLs/SortQTLtoolsPhenotypeTable/{tissue_id}.{MultiOrSingle}.{K}.log"
    resources:
        mem_mb = 24000
    wildcard_constraints:
        tissue_id = '|'.join(tissue_list),
        MultiOrSingle = 'multi_tissue|single_tissue',
        K = '2|3|4|5|10'
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
        'QTLs/{tissue_id}/{MultiOrSingle}.snmf_{K}.sorted.qqnorm.bed.gz',
    output:
        'QTLs/{tissue_id}/{MultiOrSingle}.snmf_{K}.sorted.qqnorm.bed.pca',
    log:
        "logs/QTLs/PhenotypePCs/{tissue_id}.{MultiOrSingle}.snmf_{K}.log"
    resources:
        mem_mb = 24000
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/PermuteAndPCA.R {input} {output} &> {log}
        """

rule collect_PCA:
    input:
        expand('QTLs/{tissue_id}/multi_tissue.snmf_{K}.sorted.qqnorm.bed.pca', tissue_id = tissue_sub_list, K = [2, 3, 4, 5, 10]),
        expand('QTLs/{tissue_id}/single_tissue.snmf_3.sorted.qqnorm.bed.pca', 
               tissue_id = ['Brain_Frontal_Cortex_BA9', 'Muscle_Skeletal', 'Whole_Blood', 'Skin_Not_Sun_Exposed_Suprapubic'])

N_PermutationChunks = 100
ChunkNumbers = range(0, 1+N_PermutationChunks) 


def GetQTLtoolsPassFlags(wildcards):
    if wildcards.Pass == "PermutationPass":
        return "--permute 1000"
    elif wildcards.Pass == "NominalPass":
        return "--nominal 1"

# In case I need to exclude outliers
def GetExcludeFile(wildcards):
    if wildcards.SexGroup == '.female':
        return "--exclude-samples config/female_participants.txt"
    elif wildcards.SexGroup == '.male':
        return "--exclude-samples config/male_participants.txt"
    else:
        return ""

rule QTLtools_generalized:
    input:
        vcf = "/project2/yangili1/cfbuenabadn/ChromatinSplicingQTLs/code/GTEx/data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.vcf.gz",
        tbi = "/project2/yangili1/cfbuenabadn/ChromatinSplicingQTLs/code/GTEx/data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.vcf.gz.tbi",
        bed = 'QTLs/{tissue_id}/{MultiOrSingle}.snmf_{K}.sorted.qqnorm.bed.gz',
        bed_tbi = 'QTLs/{tissue_id}/{MultiOrSingle}.snmf_{K}.sorted.qqnorm.bed.gz.tbi',
        cov = 'QTLs/{tissue_id}/{MultiOrSingle}.snmf_{K}.sorted.qqnorm.bed.pca'
    output:
        temp("QTLs/{tissue_id}/{MultiOrSingle}{SexGroup}.snmf_{K}.{Pass}.Chunks/{QTLTools_chunk_n}.txt")
    log:
        "logs/QTLs/{tissue_id}/{MultiOrSingle}{SexGroup}.snmf_{K}.{Pass}.Chunks/{QTLTools_chunk_n}.log"
    resources:
        mem_mb = much_more_mem_after_first_attempt
    envmodules:
        "gsl/2.5"
    params:
        WindowFlag = "--window 10000",
        OtherFlags = "--grp-best",
        PassFlags = GetQTLtoolsPassFlags,
        ExcFlag = GetExcludeFile,
        NChunks = N_PermutationChunks
    wildcard_constraints:
        n = "|".join(str(i) for i in ChunkNumbers),
        Pass = "PermutationPass|NominalPass",
        SexGroup = "|.female|.male",
        MultiOrSingle = 'multi_tissue|single_tissue'
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
        expand( "QTLs/{{tissue_id}}/{{MultiOrSingle}}{{SexGroup}}.snmf_{{K}}.{{Pass}}.Chunks/{QTLTools_chunk_n}.txt", QTLTools_chunk_n=ChunkNumbers )
    output:
        "QTLs/{tissue_id}/{MultiOrSingle}{SexGroup}.snmf_{K}.{Pass}.txt.gz"
    log:
        "logs/Gather_QTLtools_cis_pass/{tissue_id}.{MultiOrSingle}{SexGroup}.snmf_{K}.{Pass}.log"
    wildcard_constraints:
        Pass = "PermutationPass|NominalPass",
        SexGroup = "|.female|.male",
        MultiOrSingle = 'multi_tissue|single_tissue'
    shell:
        """
        (cat {input} | gzip - > {output}) &> {log}
        """
        

rule AddQValueToPermutationPass:
    input:
        "QTLs/{tissue_id}/{MultiOrSingle}{SexGroup}.snmf_{K}.PermutationPass.txt.gz"
    output:
        table = "QTLs/{tissue_id}/{MultiOrSingle}{SexGroup}.snmf_{K}.PermutationPass.FDR_Added.txt.gz",
    log:
        "logs/AddQValueToPermutationPass/{tissue_id}/{MultiOrSingle}{SexGroup}.snmf_{K}.PermutationPass.txt.gz.log"
    conda:
        "../envs/r_essentials.yml"
    wildcard_constraints:
        SexGroup = "|.female|.male",
        MultiOrSingle = 'multi_tissue|single_tissue'
    priority:
        10
    shell:
        """
        Rscript scripts/AddQvalueToQTLtoolsOutput.R {input} {output.table} PermutationPass &> {log}
        """

rule CollectQTLs:
    input:
        expand("QTLs/{tissue_id}/multi_tissue{SexGroup}.snmf_{K}.{Pass}.txt.gz", 
               Pass = ['PermutationPass.FDR_Added', 'NominalPass'],
               tissue_id = ['Brain_Frontal_Cortex_BA9', 'Muscle_Skeletal', 'Whole_Blood', 'Skin_Not_Sun_Exposed_Suprapubic'],
               SexGroup = ['', '.female', '.male'],
               K = [2, 3, 4, 5, 10])





rule collect_gao_EL:
    output:
        expand('ebpmf_models/gao_models/tables/snmf{k}_EL.tab.gz', k=[3, 5])
    log:
         'logs/collect_gao_EL.log'
    resources:
        mem_mb = 12000
    shell:
        """
        python scripts/collect_ad_EL.py &> {log}
        """

rule collect_gregor_EL:
    output:
        expand('ebpmf_models/sberger_models/tables/{seq}_snmf{k}_EL.tab.gz', k=[3, 5], seq=['seqmatic', 'invitae'])
    log:
         'logs/collect_sberger_EL.log'
    resources:
        mem_mb = 12000
    shell:
        """
        python scripts/collect_sberger_EL.py &> {log}
        """

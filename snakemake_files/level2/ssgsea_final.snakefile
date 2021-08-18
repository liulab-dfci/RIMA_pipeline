#!/usr/bin/env python

#-------------------------------ssGSEA-----------------------#
data_type = config["assembly"]

def DiffReference(data_type):
    if data_type == "mm10":
        reference = "static/ssgsea/kegg.v5.2.symbols_mouse.gmt"
    if data_type == "hg38":
        reference = "static/ssgsea/c2.cp.kegg.v6.1.symbols.gmt"
    return reference

def ssgsea_targets(wildcards):
    ls = []
    ls.append("files/ssgsea/ssgsea_all_terms_score.txt")
    ls.append("files/multiqc/ssgsea/ssgsea.txt")
    for comparison in config['designs']:
        ls.append("images/ssgsea/%s/%s_ssgsea_plot.png" % (comparison, comparison))
    return ls

rule ssgsea_all:
    input:
        ssgsea_targets

rule ssgsea_score:
    input:
        "analysis/combat/tpm_convertID.batch"
    output:
        score = "files/ssgsea/ssgsea_all_terms_score.txt",
        plots = expand("images/ssgsea/{comparison}/{comparison}_ssgsea_plot.png", comparison=config['designs']),
        ssgsea_heatmap = "files/multiqc/ssgsea/ssgsea.txt"
    log:
        "logs/ssgsea/ssgsea.log"
    params:
        gmt = DiffReference(data_type),
        comparison = lambda wildcards: ','.join(str(i) for i in list(config['designs'])),
        meta = config['metasheet'],
        top_n = "50",
        outpath = "images/ssgsea/",
        path="set +eu;source activate %s" % config['stat_root'],
        order= config['designs'],
        multiqc_folder = "files/multiqc/ssgsea/"
    message:
        "Running single sample gene set enrichment analysis" 
    benchmark:
        "benchmarks/ssgsea/ssgsea.benchmark"
    conda: "../envs/stat_perl_r.yml"
    shell:
        "{params.path}; Rscript src/ssgsea_plot.R -e {input} -g {params.gmt} -c {params.comparison} -m {params.meta} -n {params.top_n} -o {params.outpath} --order {params.order} --multiqc {params.multiqc_folder}"
        " && mv images/ssgsea/ssgsea_all_terms_score.txt files/ssgsea/ssgsea_all_terms_score.txt"
        

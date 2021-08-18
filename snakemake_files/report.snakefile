#!/usr/bin/env python

#-------------------------------multiqc report----------------------#
# @author: Lin Yang, Aashna
# @email:

def report_targets(wildcards):
    ls = []
    ls.append("files/multiqc/QC_check/QC_status.txt")
    ls.append("files/multiqc/QC_check/meta_information.csv")
    ls.append("files/RIMA.analysis.html")
    return ls


rule multiqc_all:
    input:
         report_targets

rule QC_Check:
    input:
        distr = "analysis/rseqc/read_distrib/read_distrib.matrix.tab",
        star = "files/star/STAR_Align_Report.csv",
        tin = "files/rseqc/tin_score/tin_score_summary.txt"
    output:
        QC_check = "files/multiqc/QC_check/QC_status.txt",
        meta_info = "files/multiqc/QC_check/meta_information.csv"
    message:
        "Checking the qulity of data"
    benchmark:
        "benchmarks/multiqc_QC/QC.benchmark"
    params:
        meta_file = config["metasheet"],
        prefix = "files/multiqc/QC_check/",
        path="set +eu;source activate %s" % config['stat_root'],
    log:
        "logs/multiqc/QC_check.log"
    conda: "../envs/stat_perl_r.yml"
    shell:
        "{params.path}; Rscript src/multiqc/QC_check.R -m {params.meta_file} -s {input.star} -d {input.distr} -t {input\
.tin} -o {params.prefix}"


rule run_multiqc:
    input:
       file="files/multiqc/"
    output:
       report= "files/RIMA.analysis.html"
    params:
       config_file= "static/multiqc/multiqc_config.yaml"
    log:
       "logs/multiqc/multiqcrun_check.log"
    shell:
       "multiqc -c {params.config_file} -n {output.report} {input.file}"

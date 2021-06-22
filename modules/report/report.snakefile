#!/usr/bin/env python

#-------------------------------multiqc report generation module----------------------#

def report_targets(wildcards):
    ls = []
    ls.append("files/multiqc/QC_check/QC_status.txt")
    ls.append("files/multiqc/QC_check/meta_information.csv")
    ls.append("files/multiqc/immune_response/tide_table.txt")
    ls.append("files/report/RIMA.analysis.html")
    return ls


rule multiqc_all:
    input:
         report_targets

rule QC_Check:
    input:
        distr = "analysis/rseqc/read_distrib/read_distrib.matrix.tab",
        star = "analysis/star/STAR_Align_Report.csv",
        tin = "analysis/rseqc/tin_score/tin_score_summary.txt"
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
        "{params.path}; Rscript src/multiqc/QC_check.R -m {params.meta_file} -s {input.star} -d {input.distr} -t {input.tin} -o {params.prefix}"

        

rule marker_comparison:
    input:
       tide= "analysis/tide/",
       msi = "files/multiqc/immune_response/msi_score.txt"
    output:
       tide_table= "files/multiqc/immune_response/tide_table.txt"
    params:
       phen = config["msi_clinical_phenotypes"],
       meta_file = config["metasheet"],
       outdir = "files/multiqc/immune_response/",
       path="set +eu;source activate %s" % config['stat_root'],
    log:
       "logs/multiqc/marker_check.log"
    shell:
        "{params.path}; Rscript src/multiqc/maker_com.R -m {params.meta_file} -s {input.msi} -t {input.tide} -p {params.phen} -o {params.outdir}"


rule run_multiqc:
    input:
       file= "files/multiqc/"
    output:
       report= "files/report/RIMA.analysis.html"
    params:
       config_file= "static/multiqc/config.yaml",
       sample_names= "static/multiqc/Multiqc_Sample_names.txt"
    log:
       "logs/multiqc/multiqcrun_check.log"
    shell:
       "multiqc -c {params.config_file} --sample-names {params.sample_names} -n {output.report} {input.file}"


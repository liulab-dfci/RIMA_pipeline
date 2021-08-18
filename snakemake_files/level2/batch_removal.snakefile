#!/usr/bin/env python

#-------------------------------limma for batch removal----------------------#
# @author: Jin Wang
# @email: jwang0611@gmail.com

import itertools

data_type = config["assembly"]

def DiffReference(data_type):
    if data_type == "mm10":
        reference = "$8"
    if data_type == "hg38":
        reference = "$10"
    return reference

def combat_targets(wildcards):
    ls = []
    ls.append("analysis/combat/tpm_convertID.txt")
    ls.append("analysis/combat/tpm_convertID.batch")
    ls.append("analysis/combat/gencode_tx2_ensemble_gene_symbol.csv")
    ls.append("images/combat/pca_plot_before.png")
#    ls.append("images/combat/pca_plot_after.png")
    return ls


rule combat_all:
    input:
         combat_targets

rule id_convert:
    input:
        "analysis/salmon/salmon_tpm_matrix.csv"
    output:
        gene_convertID = "analysis/combat/tpm_convertID.txt",
        annotation = "analysis/combat/gencode_tx2_ensemble_gene_symbol.csv"
    message:
        "Convert Ensemble ID or transcript ID to Gene Symbol from salmon matrix"
    log:
        "logs/combat/convertID.log"
    benchmark:
        "benchmarks/combat/convertID.benchmark"
    params:
        gene_prefix = "analysis/combat/tpm_convertID.txt",
        tool_type = "salmon",
        sp = DiffReference(data_type)
    shell:
        """zless {config[gtf_path]} | grep -v "#" | awk '$3=="transcript"' | cut -f9 | tr -s ";" " \
" | awk '{{print$4"\t"$2"\t"{params.sp}}}' | sort | uniq | sed 's/\"//g' > {output.annotation} """
        """ && python src/ConvertID.py -m {input} -t {params.tool_type} -r {output.annotation} -p {\
params.gene_prefix} """


rule batch_removal:
    input:
        "analysis/combat/tpm_convertID.txt"
    output:
        combat_expr = "analysis/combat/tpm_convertID.batch"
    message:
        "Running batch removal using limma method"
    benchmark:
        "benchmarks/combat/tpm_combat.benchmark"
    params:
        batch_file = config["metasheet"],
        covariates = lambda wildcards: ','.join(str(i) for i in config["batch_covariates"]),
        prefix = "analysis/combat/tpm_convertID",
        path="set +eu;source activate %s" % config['stat_root']
    log:
        "logs/combat/batch_removal.log"
    conda: "../envs/stat_perl_r.yml"
    shell:
        "{params.path}; Rscript src/batch_removal.R -e {input} -b {params.batch_file} -c {params.covariates} -o {params.prefix}"

rule pca_sample_clustering:
    input:
        before_batch = "analysis/combat/tpm_convertID.txt",
        after_batch = "analysis/combat/tpm_convertID.batch"
    output:
        "images/combat/pca_plot_before.png",
 #       "images/combat/pca_plot_after.png"
    message:
        "Running PCA for sample clustering"
    benchmark:
        "benchmarks/combat/pca.benchmark"
    params:
        meta_info = config["metasheet"],
        out_path = "images/combat/",
        path="set +eu;source activate %s" % config['stat_root'],
        covariates = lambda wildcards: ','.join(str(i) for i in config["batch_covariates"]),
        design = config["designs"]
#        batch = config["batch_name"]
    conda: "../envs/stat_perl_r.yml"
    shell:
        "{params.path}; Rscript src/pca.R -b {input.before_batch} -a {input.after_batch} -m {params\
.meta_info} -o {params.out_path} -c {params.covariates} --design {params.design}"

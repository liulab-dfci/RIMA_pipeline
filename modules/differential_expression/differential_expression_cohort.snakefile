#!/usr/bin/env python

#-------------------------------differential gene expression-----------------------#

import itertools
import pandas as pd

metadata = pd.read_csv(config["metasheet"], index_col=0, sep=',')
options = [config["Treatment"],config["Control"]]
design = config["design"]
treatment = config["Treatment"]
control = config["Control"] 

data_type = config["assembly"]

#def DiffInput(data_type):
#        return expand("analysis/salmon/{sample}/{sample}.quant.sf",sample=config["samples"])

def DiffReference(data_type):
    if data_type == "mm10":
        reference = ["static/deseq2/tx2gene_mouse.csv", "static/deseq2/mm10_pcoding.csv"]
    if data_type == "hg38":
        reference = ["static/deseq2/tx2gene.csv", "static/deseq2/hg_pcoding.csv"]
    return reference

def ssgseaRef(data_type):
    if data_type == "mm10":
        reference = "static/ssgsea/kegg.v5.2.symbols_mouse.gmt"
    if data_type == "hg38":
        reference = "static/ssgsea/c2.cp.kegg.v6.1.symbols.gmt"
    return reference


def diffexpr_targets(wildcards):
    ls=[]
    ls.append("analysis/deseq2/%s_%s_vs_%s_volcano_plot.png" % (design,treatment,control))
    ls.append("analysis/deseq2/%s_%s_vs_%s_DESeq2.txt" % (design,treatment,control))
    ls.append("analysis/deseq2/%s_%s_vs_%s_TPMs.txt" % (design,treatment,control))
    ls.append("analysis/gsea/%s_%s_vs_%s_BP_terms.txt" % (design,treatment,control))
    ls.append("analysis/gsea/%s_%s_vs_%s_MF_terms.txt" % (design,treatment,control))
    ls.append("analysis/gsea/%s_%s_vs_%s_CC_terms.txt" % (design,treatment,control))
#    ls.append("analysis/gsea/%s_%s_vs_%s_KEGG_terms.txt" % (design,treatment,control))
    ls.append("analysis/gsea/%s_%s_vs_%s_BP_terms.png" % (design,treatment,control))
    ls.append("analysis/gsea/%s_%s_vs_%s_MF_terms.png" % (design,treatment,control))
    ls.append("analysis/gsea/%s_%s_vs_%s_CC_terms.png" % (design,treatment,control))
#    ls.append("analysis/gsea/%s_%s_vs_%s_KEGG_terms.png" % (design,treatment,control))
#    ls.append("analysis/ssgsea/%s_%s_vs_%s_kegg_ssgsea.txt" % (design,treatment,control))

    return ls

def getsampleIDs(meta):
	return meta[meta[design].isin(options)].index


rule differential_all:
    input:
        diffexpr_targets

rule deseq2_differential_genes:
    input:
    	files = expand("analysis/salmon/{sample}/{sample}.quant.sf",sample=getsampleIDs(metadata))
    output:
        deseq2_res = "analysis/deseq2/{design}_{treatment}_vs_{control}_DESeq2.txt",
        deseq2_tpm = "analysis/deseq2/{design}_{treatment}_vs_{control}_TPMs.txt"
    log:
        "logs/deseq2/{design}.{treatment}.{control}.deseq2.log"
    params:
    	filelist = lambda wildcards, input: ','.join(str(i) for i in list({input.files})[0]),
        batch = config["batch"],
        out_path = "analysis/deseq2/",
        tx_annot = DiffReference(data_type)[0],
        condition = design,
        treatment = treatment,
        control = control,
        meta = config["metasheet"],
        pgenes = DiffReference(data_type)[1],
        path = "set +eu;source activate %s" % config['stat_root']
    message:
        "Running DESeq2 on the samples"
    benchmark:
        "benchmarks/deseq2/{design}.{treatment}.{control}.deseq2.benchmark"
    conda: "../envs/stat_perl_r.yml"
    shell:
    	"{params.path}; Rscript src/differentialexpr/DESeq2.R --input {params.filelist} --type salmon \
        --batch {params.batch} --meta {params.meta} --tx2gene {params.tx_annot} \
        --condition {params.condition} --pcoding {params.pgenes} --treatment {params.treatment} --control {params.control} --outpath {params.out_path}"
        
        


rule volcano_plot:
    input:
        "analysis/deseq2/{design}_{treatment}_vs_{control}_DESeq2.txt"
    output:
        "analysis/deseq2/{design}_{treatment}_vs_{control}_volcano_plot.png"
    params:
        path = "set +eu;source activate %s" % config['stat_root'],
        treatment = treatment,
        control = control,
        condition = design
    log:
        "logs/deseq2/{design}_{treatment}_vs_{control}.deseq2.volcano.log"
    message:
        "Running Volcano plot on the DESeq2 result"
    benchmark:
        "benchmarks/deseq2/{design}_{treatment}_vs_{control}.deseq2.volcano.benchmark"
    conda: "../envs/stat_perl_r.yml"
    shell:
        "{params.path}; Rscript src/differentialexpr/volcano_plot.R --deseq2_mat {input} \
        --treatment {params.treatment} --control {params.control} --condition {params.condition} \
        --output {output}"




rule gsea_plot:
    input: "analysis/deseq2/{design}_{treatment}_vs_{control}_DESeq2.txt"
    output:
        go_bp = "analysis/gsea/{design}_{treatment}_vs_{control}_BP_terms.txt",
        go_mf = "analysis/gsea/{design}_{treatment}_vs_{control}_MF_terms.txt",
        go_cc = "analysis/gsea/{design}_{treatment}_vs_{control}_CC_terms.txt",
       # kegg = "analysis/gsea/{design}_{treatment}_vs_{control}_KEGG_terms.txt",
        go_bp_plot = "analysis/gsea/{design}_{treatment}_vs_{control}_BP_terms.png",
        go_mf_plot = "analysis/gsea/{design}_{treatment}_vs_{control}_MF_terms.png",
        go_cc_plot = "analysis/gsea/{design}_{treatment}_vs_{control}_CC_terms.png"
 #       kegg_plot = "analysis/gsea/{design}_{treatment}_vs_{control}_KEGG_terms.png",
    log:
        "logs/gsea/{design}_{treatment}_vs_{control}.gsea.log"
    params:
        out_path = "analysis/gsea/",
        gsea_pcut = 10,
        gsea_minsize = 3,
        gsea_permutation = 1000,
        treatment = treatment,
        control = control,
        condition = design,
        path = "set +eu;source activate %s" % config['stat_root']
    message:
        "Running GSEA on the samples"
    benchmark:
        "benchmarks/gsea/{design}_{treatment}_vs_{control}.gsea.benchmark"
    conda: "../envs/stat_perl_r.yml"
    shell:
        "{params.path}; Rscript src/differentialexpr/gsea.R --deseq2_mat {input} \
        --pcut {params.gsea_pcut} --minsize {params.gsea_minsize} --npermutation {params.gsea_permutation} \
        --outdir {params.out_path} --treatment {params.treatment} --control {params.control} \
        --condition {params.condition}"



#rule ssgsea:
#    input:
#        "analysis/deseq2/{design}_{treatment}_vs_{control}_TPMs.txt"
#    output:
#        score = "analysis/ssgsea/{design}_{treatment}_vs_{control}_kegg_ssgsea.txt"
#    log:
#        "logs/ssgsea/{design}_{treatment}_vs_{control}_ssgsea.log"
#    params:
#        gmt = ssgseaRef(data_type),
#        condition = design,
#        treatment = treatment,
#        control = control,
#        meta = config['metasheet'],
#        top_n = 20,
#        outpath = "analysis/ssgsea/",
#        path = "set +eu;source activate %s" % config['stat_root']
#    message:
#        "Running single sample gene set enrichment analysis"
#    benchmark:
#        "benchmarks/ssgsea/{design}_{treatment}_vs_{control}_ssgsea.benchmark"
#    conda: "../envs/stat_perl_r.yml"
#    shell:
#        "{params.path}; Rscript src/differentialexpr/ssgsea.R -e {input} -f {params.gmt} \
#        --treatment {params.treatment} --control {params.control} --condition {params.condition}\
#        -m {params.meta} -n {params.top_n} -o {params.outpath}"





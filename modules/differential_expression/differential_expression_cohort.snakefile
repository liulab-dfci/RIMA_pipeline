#!/usr/bin/env python

#-------------------------------differential gene expression-----------------------#

import itertools
import pandas as pd


def diffexpr_targets(wildcards):
    ls=[]
    ls.append("files/multiqc/differentialexpr/ssgsea.txt")
    #ls.append("files/multiqc/differentialexpr/deseq2_description.txt")
    for design in config["designs"]:
        meta = config["metasheet"]
        design_file = "./" + meta
        design_meta = pd.read_csv(design_file, index_col=0, sep=',')
        design_meta = design_meta.rename(columns = {design: 'Condition'})
        comps = [i for i in list(set(design_meta["Condition"].tolist()))]
        combinations = list(itertools.combinations(comps,2))
        if len(comps) >1:
            if config["comparison"] == "between":
                compare_list = combinations
            if config["comparison"] == "loop":
                compare_list = comps
            for cp in compare_list:
                if config["comparison"] == "between":
                    cp = sorted(list(cp))
                    cp_list = cp[0]+"_VS_"+cp[1]
                else:
                    cp_list = cp+"_VS_others"
                ls.append("files/deseq2/%s/%s_%s_diff_volcano_plot.png" % (design,design,cp_list))
                ls.append("analysis/deseq2/%s/%s_%s_DESeq2_ConvertID.txt" % (design,design,cp_list))
                ls.append("files/multiqc/differentialexpr/%s/%s_%s_GO_BP_terms.txt" % (design,design,cp_list))
                ls.append("files/multiqc/differentialexpr/%s/%s_%s_GO_CC_terms.txt" % (design,design,cp_list))
                ls.append("files/multiqc/differentialexpr/%s/%s_%s_GO_MF_terms.txt" % (design,design,cp_list))
                ls.append("files/multiqc/differentialexpr/%s/%s_%s_KEGG_terms.txt" % (design,design,cp_list))
                ls.append("files/multiqc/differentialexpr/%s_%s_DESeq2_sub.txt" % (design,cp_list))
    return ls

rule differential_all:
    input:
        diffexpr_targets

rule deseq2_differential_genes:
    input:
        files = expand("analysis/salmon/{sample}/{sample}.quant.sf",sample=config["samples"]),
    output:
        deseq2_res = "analysis/deseq2/{design}/{design}_{compare}_DESeq2_ConvertID.txt",
        multiqc = "files/multiqc/differentialexpr/{design}_{compare}_DESeq2_sub.txt"
    log:
        "logs/deseq2/{design}/{design}_{compare}.deseq2.log"
    params:
        comps = config["comparison"],
        filelist = lambda wildcards, input: ','.join(str(i) for i in list({input.files})[0]),
        batch = config["batch"],
        out_path = "analysis/deseq2/{design}/{design}",
        file_path = "files/deseq2/{design}/",
        multiqc_path = "files/multiqc/differentialexpr/deseq2_description.txt",
        tx_annot = "static/deseq2/tx2gene.csv",
        condition = config["designs"],
        meta = config["metasheet"],
        source = "salmon",
        path = "set +eu;source activate %s" % config['stat_root']
    message:
        "Running DESeq2 on the samples"
    benchmark:
        "benchmarks/deseq2/{design}/{design}_{compare}.deseq2.benchmark"
    conda: "../envs/stat_perl_r.yml"
    shell:
        "{params.path}; Rscript src/differentialexpr/DESeq2.R --input {params.filelist} --batch {params.batch} --type {params.source} --comparison {params.comps} --meta {params.meta} --tx2gene {params.tx_annot} -\
-condition {params.condition} --outpath {params.out_path} --multiqc {output.multiqc}"
        " && cp {output.multiqc} {params.multiqc_path}"
        #" && cp {params.out_path}*_DESeq2_ConvertID.txt {params.file_path}"
        #" && cp {params.out_path}*_raw.txt {params.multiqc_path}"


rule volcano_plot:
    input:
        "analysis/deseq2/{design}/{design}_{compare}_DESeq2_ConvertID.txt"
    output:
        "files/deseq2/{design}/{design}_{compare}_diff_volcano_plot.png"
    params:
        path = "set +eu;source activate %s" % config['stat_root'],
        multiqc_path = "files/multiqc/differentialexpr/DESeq2-Volcano-Plot_mqc.png",
        meta = config["metasheet"],
        condition = config["designs"]
    log:
        "logs/deseq2/{design}/{design}_{compare}.deseq2.volcano.log"
    message:
        "Running Volcano plot on the DESeq2 result"
    benchmark:
        "benchmarks/deseq2/{design}/{design}_{compare}.deseq2.volcano.benchmark"
    conda: "../envs/stat_perl_r.yml"
    shell:
        "{params.path}; Rscript src/differentialexpr/volcano_plot.R --deseq2_mat {input} --meta {params.meta} --condition {params.condition} --outdir {output}"
        " && cp {output} {params.multiqc_path} "


rule gsea_plot:
    input: "analysis/deseq2/{design}/{design}_{compare}_DESeq2_ConvertID.txt"
    output:
        go_bp = "files/multiqc/differentialexpr/{design}/{design}_{compare}_GO_BP_terms.txt",
        go_mf = "files/multiqc/differentialexpr/{design}/{design}_{compare}_GO_MF_terms.txt",
        go_cc = "files/multiqc/differentialexpr/{design}/{design}_{compare}_GO_CC_terms.txt",
        kegg = "files/multiqc/differentialexpr/{design}/{design}_{compare}_KEGG_terms.txt"
    log:
        "logs/gsea/{design}/{design}_{compare}.gsea.log"
    params:
        out_path = "files/multiqc/differentialexpr/{design}/{design}_{compare}",
        gsea_pcut = "0.1",
        gsea_minsize = "5",
        gsea_permutation = "1000",
        path = "set +eu;source activate %s" % config['stat_root']
    message:
        "Running GSEA on the samples"
    benchmark:
        "benchmarks/gsea/{design}/{design}_{compare}.gsea.benchmark"
    conda: "../envs/stat_perl_r.yml"
    shell:
        "{params.path}; Rscript src/differentialexpr/gsea.R --deseq2_mat {input} --pcut {params.gsea_pcut} --minsize {params.gsea_minsize} --npermutation {params.gsea_permutation} -o {params.out_path}"
        #" && mv {params.out_path}*txt {params.files_path}"

rule ssgsea_score:
    input:
        "analysis/batchremoval/tpm_convertID.batch"
    output:
        score = "files/multiqc/differentialexpr/ssgsea.txt",
    log:
        "logs/ssgsea/ssgsea.log"
    params:
        gmt = config['gmt_path'],
        comparison = lambda wildcards: ','.join(str(i) for i in list(config['ssgsea_comparisons'])),
        meta = config['metasheet'],
        top_n = 50,
        outpath = "files/multiqc/differentialexpr/",
        path = "set +eu;source activate %s" % config['stat_root'],
        order = config['sam_ord']
    message:
        "Running single sample gene set enrichment analysis"
    benchmark:
        "benchmarks/ssgsea/ssgsea.benchmark"
    conda: "../envs/stat_perl_r.yml"
    shell:
        "{params.path}; Rscript src/differentialexpr/ssgsea.R -e {input} -g {params.gmt} -c {params.comparison} -m {params.meta} -n {params.top_n} -o {params.outpath} --order {params.order}"

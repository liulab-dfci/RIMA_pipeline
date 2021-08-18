#!/usr/bin/env python

#-------------------------------GSEA for pathway enrichment analysis-----------------------#

import itertools
import pandas as pd

def gsea_targets(wildcards):
    ls=[]
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
                ls.append("images/gsea/%s/%s_%s_diff_gsea_plot.png" % (design,design,cp_list))
                ls.append("files/gsea/%s/%s_%s_GO_BP_terms.txt" % (design,design,cp_list))
                ls.append("files/gsea/%s/%s_%s_GO_BP_terms.txt" % (design,design,cp_list))
                ls.append("files/gsea/%s/%s_%s_GO_BP_terms.txt" % (design,design,cp_list))
                ls.append("files/gsea/%s/%s_%s_KEGG_terms.txt" % (design,design,cp_list))
    return ls



rule gsea_plot:
    input: "files/deseq2/{design}/{design}_{compare}_DESeq2_ConvertID.txt"
    output:
        gsea_plot = "images/gsea/{design}/{design}_{compare}_diff_gsea_plot.png",
        go_bp = "files/gsea/{design}/{design}_{compare}_GO_BP_terms.txt",
        go_mf = "files/gsea/{design}/{design}_{compare}_GO_MF_terms.txt",
        go_cc = "files/gsea/{design}/{design}_{compare}_GO_CC_terms.txt",
        kegg = "files/gsea/{design}/{design}_{compare}_KEGG_terms.txt"
    log:
        "logs/gsea/{design}/{design}_{compare}.gsea.log"
    params:
        files_path = "files/gsea/{design}/",
        out_path = "images/gsea/{design}/{design}_{compare}",
        gsea_pcut = "0.1",
        gsea_species = config['assembly'],
        gsea_permutation = "1000",
        path="set +eu;source activate %s" % config['stat_root']
    message:
        "Running GSEA on the samples"
    benchmark:
        "benchmarks/gsea/{design}/{design}_{compare}.gsea.benchmark"
    conda: "../envs/stat_perl_r.yml"
    shell:
        "{params.path}; Rscript src/gsea_plot.R --deseq2_mat {input} --pcut {params.gsea_pcut} --sp\
ecies {params.gsea_species} --npermutation {params.gsea_permutation} --outdir {params.out_path}"
        " && mv {params.out_path}*txt {params.files_path}"


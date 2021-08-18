#!/usr/bin/env python

#-------------------------------
# @author: Jin Wang; @ref: VIPER
# @email: jwang0611@gmail.com


def immunedeconv_targets(wildcards):
    ls = []
    ls.append("files/immunedeconv/quantiseq.txt")
    ls.append("files/immunedeconv/xcell.txt")
    ls.append("files/immunedeconv/mcp_counter.txt")
    ls.append("files/immunedeconv/cibersort_abs.txt")
    ls.append("files/immunedeconv/timer.txt")
    ls.append("files/immunedeconv/epic.txt")
    ls.append("images/immunedeconv/B_cell_corr.png")
    ls.append("images/immunedeconv/DC_corr.png")
    ls.append("images/immunedeconv/NK_corr.png")
    ls.append("images/immunedeconv/CD4_T_cell_corr.png")
    ls.append("images/immunedeconv/Neutrophil_corr.png")
    ls.append("images/immunedeconv/CD8_T_cell_corr.png")
    ls.append("images/immunedeconv/Macrophage_corr.png")
    ls.append("images/immunedeconv/Treg_corr.png")
    ls.append("images/immunedeconv/ImmuneDeconv_heatmpap.png")
    ls.append("files/multiqc/immunedeconv/cibersort_abs_immune_multiqc.txt")
    return ls

rule immunedeconv_all:
    input:
        immunedeconv_targets

rule ImmuneDeconv_infiltration:
    input:
        "analysis/combat/tpm_convertID.txt"
    output:
        "files/immunedeconv/quantiseq.txt",
        "files/immunedeconv/xcell.txt",
        "files/immunedeconv/mcp_counter.txt",
        "files/immunedeconv/cibersort_abs.txt",        
        "files/immunedeconv/timer.txt",
        "files/immunedeconv/epic.txt"
    log:
        "logs/immunedeconv/deconv.log"
    params:
        perm = "100",
        qn = "FALSE",
        absl = "TRUE",
        abs_method = "sig.score",
        out_dir = "files/immunedeconv",
        ref = "static/immunedeconv/mm_to_hg.csv",
        path="set +eu;source activate %s" % config['immnue_decov']
    message: 
        "Running immunedeconv on the expression data"
    benchmark:
        "benchmarks/immunedeconv/deconv.benchmark"
    shell:
        "{params.path}; Rscript src/immunedeconv.R -e {input} -t {config[cancer_type]} -p {params.perm} -q {params.qn} -a {params.absl} -m {params.abs_method} -o {params.out_dir} -d {config[assembly]} -r {params.ref}"

rule ImmuneDeconv_plot:
    input:
        "files/immunedeconv/quantiseq.txt",
        "files/immunedeconv/xcell.txt",
        "files/immunedeconv/mcp_counter.txt",
        "files/immunedeconv/cibersort_abs.txt",        
        "files/immunedeconv/timer.txt",
        "files/immunedeconv/epic.txt"
    output:
        "images/immunedeconv/B_cell_corr.png",
        "images/immunedeconv/DC_corr.png",
        "images/immunedeconv/NK_corr.png",
        "images/immunedeconv/CD4_T_cell_corr.png",
        "images/immunedeconv/Neutrophil_corr.png",
        "images/immunedeconv/CD8_T_cell_corr.png",
        "images/immunedeconv/Macrophage_corr.png",
        "images/immunedeconv/Treg_corr.png",
        "images/immunedeconv/ImmuneDeconv_heatmpap.png",
        "files/multiqc/immunedeconv/cibersort_abs_immune_multiqc.txt"
    log:
        "logs/immunedeconv/deconv_plot.log"
    message: 
        "Running immunedeconv plot"
    benchmark:
        "benchmarks/immunedeconv/deconv_plot.benchmark"
    conda: "../envs/stat_perl_r.yml"
    params:
        meta = config["metasheet"],
        inputdir = "files/immunedeconv/",
        outdir = "images/immunedeconv/",
        multiqc = "files/multiqc/immunedeconv/",
        clinic_phenotype = config["designs"],
        path="set +eu;source activate %s" % config['stat_root']
    shell:
        "{params.path}; Rscript src/immunedeconv_plot.R --meta {params.meta} --input_dir {params.inputdir} --output_dir {params.outdir} --clinic_col {params.clinic_phenotype} --multiqc {params.multiqc}"
    




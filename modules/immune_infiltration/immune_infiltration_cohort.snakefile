#!/usr/bin/env python

#metadata = pd.read_csv(config["metasheet"], index_col=0, sep=',')
#options = [config["Treatment"],config["Control"]]
design = config["design"]
#treatment = config["Treatment"]
#control = config["Control"] 
batch = config["batch"]

############################Immune Infiltration by Immunedeconv######################3
def immune_infiltration_targets(wildcards):
    ls = []
    ls.append("analysis/immune_infiltration/%s_%s_quantiseq.txt" % (design,batch))
    ls.append("analysis/immune_infiltration/%s_%s_xcell.txt" % (design,batch))
    ls.append("analysis/immune_infiltration/%s_%s_mcp_counter.txt" % (design,batch))
    ls.append("analysis/immune_infiltration/%s_%s_cibersort_abs.txt" % (design,batch))
    ls.append("analysis/immune_infiltration/%s_%s_timer.txt" % (design,batch))
    ls.append("analysis/immune_infiltration/%s_%s_epic.txt" % (design,batch))
    ls.append("analysis/immune_infiltration/%s_%s_B_cell_corr.png" % (design,batch))
    ls.append("analysis/immune_infiltration/%s_%s_DC_corr.png" % (design,batch))
    ls.append("analysis/immune_infiltration/%s_%s_NK_corr.png" % (design,batch))
    ls.append("analysis/immune_infiltration/%s_%s_CD4_T_cell_corr.png" % (design,batch))
    ls.append("analysis/immune_infiltration/%s_%s_Neutrophil_corr.png" % (design,batch))
    ls.append("analysis/immune_infiltration/%s_%s_CD8_T_cell_corr.png" % (design,batch))
    ls.append("analysis/immune_infiltration/%s_%s_Macrophage_corr.png" % (design,batch))
    ls.append("analysis/immune_infiltration/%s_%s_Treg_corr.png" % (design,batch))
    #ls.append("analysis/immune_infiltration/%s_%s_ImmuneDeconv_heatmap.png" % (design,batch))
    return ls

rule immune_infiltration__all:
    input:
        immune_infiltration_targets

rule ImmuneDeconv_infiltration:
    input:
        "analysis/batchremoval/{design}_{batch}_tpm.genesymbol.batchremoved.csv" 
    output:
        "analysis/immune_infiltration/{design}_{batch}_quantiseq.txt",
        "analysis/immune_infiltration/{design}_{batch}_xcell.txt",
        "analysis/immune_infiltration/{design}_{batch}_mcp_counter.txt",
        "analysis/immune_infiltration/{design}_{batch}_cibersort_abs.txt",
        "analysis/immune_infiltration/{design}_{batch}_timer.txt",
        "analysis/immune_infiltration/{design}_{batch}_epic.txt"
    log:
        "logs/immune_infiltration/{design}_{batch}_deconv.log"
    params:
        perm = 100,
        qn = "FALSE",
        absl = "TRUE",
        abs_method = "sig.score",
        out_dir = "analysis/immune_infiltration/{design}_{batch}_",
        ref = "static/immunedeconv/mm_to_hg.csv",
        path="set +eu;source activate %s" % config['stat_root']
    message:
        "Running immunedeconv on the expression data"
    benchmark:
        "benchmarks/immune_infiltration/{design}_{batch}_deconv.benchmark"
    conda: "../envs/stat_perl_r.yml"
    shell:
        "{params.path}; Rscript src/immune_infiltration/immune_infiltration.R -e {input} -t {config[cancer_type]} -p {params.perm} -q {params.qn} \
                -a {params.absl} -m {params.abs_method} -o {params.out_dir} -d {config[assembly]} -r {params.ref}"


rule ImmuneDeconv_plot:
    input:
        "analysis/immune_infiltration/{design}_{batch}_quantiseq.txt",
        "analysis/immune_infiltration/{design}_{batch}_xcell.txt",
        "analysis/immune_infiltration/{design}_{batch}_mcp_counter.txt",
        "analysis/immune_infiltration/{design}_{batch}_cibersort_abs.txt",
        "analysis/immune_infiltration/{design}_{batch}_timer.txt",
        "analysis/immune_infiltration/{design}_{batch}_epic.txt"
    output:
        "analysis/immune_infiltration/{design}_{batch}_B_cell_corr.png",
        "analysis/immune_infiltration/{design}_{batch}_DC_corr.png",
        "analysis/immune_infiltration/{design}_{batch}_NK_corr.png",
        "analysis/immune_infiltration/{design}_{batch}_CD4_T_cell_corr.png",
        "analysis/immune_infiltration/{design}_{batch}_Neutrophil_corr.png",
        "analysis/immune_infiltration/{design}_{batch}_CD8_T_cell_corr.png",
        "analysis/immune_infiltration/{design}_{batch}_Macrophage_corr.png",
        "analysis/immune_infiltration/{design}_{batch}_Treg_corr.png",
        #"analysis/immune_infiltration/{design}_{batch}_ImmuneDeconv_heatmap.png",
    log:
        "logs/immune_infiltration/{design}_{batch}_deconv_plot.log"
    message:
        "Running immunedeconv plot"
    benchmark:
        "benchmarks/immune_infiltration/{design}_{batch}_deconv_plot.benchmark"
    conda: "../envs/stat_perl_r.yml"
    params:
        meta = config["metasheet"],
        inputdir = "analysis/immune_infiltration/{design}_{batch}_",
        outdir = "analysis/immune_infiltration/{design}_{batch}_",
        clinic_phenotype = config["design"],
        path="set +eu;source activate %s" % config['stat_root']
    shell:
        "{params.path}; Rscript src/immune_infiltration/immune_infiltration_plot.R --meta {params.meta} --input_dir {params.inputdir} --output_dir {params.outdir} --clinic_col {params.clinic_phenotype}"

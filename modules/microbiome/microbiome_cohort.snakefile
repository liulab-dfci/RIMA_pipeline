#!/usr/bin/env python

#-------------------------------Microbiota Cohort Classification---------------------------#
import pandas as pd

metadata = pd.read_csv(config["metasheet"], index_col=0, sep=',')
options = [config["Treatment"],config["Control"]]
design = config["design"]

def getsampleIDs(meta):
	return meta[meta[design].isin(options)].index
	

def microbiome_cohort_targets(wildcards):
    ls = []
    ls.append("analysis/microbiome/%s_merged_microbiota_abundance.txt" % design)
    ls.append("analysis/microbiome/%s_selected_microbes_taxonID.txt" % design)
    ls.append("analysis/microbiome/%s_Microbiome_abundance_ratio.txt" % design)
    return ls

rule microbiome_cohort_all:
    input:
      microbiome_cohort_targets

rule microbiome_plot:
    input:
      expand( "analysis/microbiome/{sample}/{sample}_addSample_report.txt", sample=getsampleIDs(metadata) ),
    output:
      merged_file = "analysis/microbiome/{design}_merged_microbiota_abundance.txt",
      selected_ID = "analysis/microbiome/{design}_selected_microbes_taxonID.txt",
      multiqc = "analysis/microbiome/{design}_Microbiome_abundance_ratio.txt"
    params:
      outdir = "analysis/microbiome/{design}_",
      phenotype = config['design'],
      meta = config["metasheet"],
      path="set +eu;source activate %s" % config['stat_root'],
    message:
      "Running Microbiota plot"
    benchmark:
      "benchmarks/centrifuge/{design}_centrifuge_merge.benchmark"
    conda: "../envs/stat_perl_r.yml"
    shell:
      "cat {input} |  sed '1 !{{/sample/d;}}' > {output.merged_file}"
      " && {params.path}; Rscript src/microbiome/mic_plot.R --input {output.merged_file} \
      --outdir {params.outdir} --clinic_col {params.phenotype} --meta {params.meta} "

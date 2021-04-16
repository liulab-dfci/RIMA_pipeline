#!/usr/bin/env python

#-------------------------------Microbiota Cohort Classification---------------------------#
import pandas as pd

metadata = pd.read_csv(config["metasheet"], index_col=0, sep=',')


def microbiome_cohort_targets(wildcards):
    ls = []
    ls.append("files/microbiome/merged_microbiota_abundance.txt")
    ls.append("files/microbiome/selected_microbes_taxonID.txt")
    ls.append("files/multiqc/microbiome/Microbiome_abundance_ratio.txt")
    return ls

rule microbiome_cohort_all:
    input:
      microbiome_cohort_targets

rule microbiome_plot:
    input:
      expand( "analysis/microbiome/{sample}/{sample}_addSample_report.txt", sample=config["samples"] ),
    output:
      merged_file = "files/microbiome/merged_microbiota_abundance.txt",
      selected_ID = "files/microbiome/selected_microbes_taxonID.txt",
      multiqc = "files/multiqc/microbiome/Microbiome_abundance_ratio.txt"
    params:
      outdir = "files/microbiome/",
      phenotype = config['sam_ord'],
      mul_folder = "files/multiqc/microbiome/",
      meta = config["metasheet"],
      path="set +eu;source activate %s" % config['stat_root'],
    message:
      "Running Microbiota plot"
    benchmark:
      "benchmarks/centrifuge/centrifuge_merge.benchmark"
    conda: "../envs/stat_perl_r.yml"
    shell:
      "cat {input} |  sed '1 !{{/sample/d;}}' > {output.merged_file}"
      " && mkdir -p {params.outdir} "
      " && {params.path}; Rscript src/microbiome/mic_plot_final.R --input {output.merged_file} --outdir {params.outdir} --clinic_col {params.phenotype} --meta {params.meta} --multiqc {params.mul_folder}"

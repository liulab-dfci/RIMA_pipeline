#!/usr/bin/env python

#-------------------------------STAR Fusion cohort-----------------------#
metadata = pd.read_csv(config["metasheet"], index_col=0, sep=',')


def mutation_cohort_targets(wildcards):
    ls = []
    ls.append("files/fusion/merged_predictions.abridged_addSample.tsv"),
    ls.append("files/fusion/pyprada_fusion_table.txt"),
    ls.append("files/fusion/pyprada_output.txt"),
    ls.append("files/multiqc/mutation/fusion/fusion_gene_table.txt"),
    ls.append("files/multiqc/mutation/fusion/fusion_gene_plot.png"),
    ls.append("files/multiqc/mutation/fusion/prada_homology.png")
    return ls

rule mutation_cohort:
    input:
      mutation_cohort_targets

rule merge_starfusion:
    input:
      expand("analysis/fusion/{sample}/{sample}.fusion_predictions.abridged_addSample.tsv", sample=config["samples"])
    output:
      "files/fusion/merged_predictions.abridged_addSample.tsv"
    log:
      "logs/fusion/merge_fusion.log"
    message:
      "Merging fusion results"
    benchmark:
      "benchmarks/fusion/merge_fusion.benchmark"
    conda: "../envs/stat_perl_r.yml"
    params:
      meta = config['metasheet']
    shell:
      """cat {input} | sed '1 !{{/FusionName/d;}}' | awk '{{gsub(/#/,"",$1);print $0}}' > {output[0]} """

rule preprocess_prada:
    input:
      "files/fusion/merged_predictions.abridged_addSample.tsv"
    output:
      "files/fusion/pyprada_fusion_table.txt"
    log:
      "logs/fusion/preprocess_prada.log"
    message:
      "Preprocessing prada input file"
    benchmark:
      "benchmarks/fusion/preprocess_prada.benchmark"
    conda: "../envs/stat_perl_r.yml"
    params:
      outdir = "files/fusion/",
      path = "set +eu;source activate %s" % config['stat_root'],
    shell:
      "{params.path}; Rscript src/mutation/preprocess_prada.R --fusion {input}  --outdir {params.outdir}"


rule run_prada:
    input:
      "files/fusion/pyprada_fusion_table.txt"
    output:
      "files/fusion/pyprada_output.txt"
    log:
      "logs/fusion/run_prada.log"
    message:
      "Running pyprada"
    benchmark:
      "benchmarks/fusion/pyprada.benchmark"
    conda: "../envs/stat_perl_r.yml"
    params:
      outdir = "files/fusion",
      config = "static/fusion/prada_config.txt",
      tmpout= "files/fusion/tmp/",
      prada_path = config["prada_path"],
      path = "set +eu;source activate %s" % config['prada_root']
    shell:
      "{params.path}"
      "&& {params.prada_path}/prada-homology -i {input} -o {output} -tmpdir {params.tmpout} -conf {params.config}"

rule fusion_plot:
    input:
      prada_input = "files/fusion/pyprada_output.txt",
      fusion = "files/fusion/merged_predictions.abridged_addSample.tsv",
      tpm_batch = "analysis/batchremoval/tpm_convertID.batch"
    output:
      fusion_table = "files/multiqc/mutation/fusion/fusion_gene_table.txt",
      fusion_plot = "files/multiqc/mutation/fusion/fusion_gene_plot.png",
      prada_plot = "files/multiqc/mutation/fusion/prada_homology.png",
    log:
      "logs/fusion/fusion_plot.log"
    message:
      "Running fusion plotting"
    benchmark:
      "benchmarks/fusion/fusion_plot.benchmark"
    conda: "../envs/stat_perl_r.yml"
    params:
      outdir = "files/multiqc/mutation/fusion/",
      prada_path=config["prada_path"],
      meta = config["metasheet"],
      path = "set +eu;source activate %s" % config['stat_root'],
      annotation = "static/fusion/cancerGeneList.tsv",
      phenotype = lambda wildcards: ','.join(str(i) for i in config["fusion_clinical_phenotypes"])
    shell:
      "{params.path}; Rscript src/mutation/fusion_plot.R  --pradafusion {input.prada_input}  --meta {params.meta}  --expression {input.tpm_batch}  --annot  {params.annotation}  --outdir  {params.outdir}  --phenotype {params.phenotype}  --input {input.fusion}"

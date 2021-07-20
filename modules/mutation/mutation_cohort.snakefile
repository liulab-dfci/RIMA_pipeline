#!/usr/bin/env python

#-------------------------------STAR Fusion cohort-----------------------#
metadata = pd.read_csv(config["metasheet"], index_col=0, sep=',')

metadata = pd.read_csv(config["metasheet"], index_col=0, sep=',')
options = [config["Treatment"],config["Control"]]
design = config["design"]
#batch = config["batch"]

def getsampleIDs(meta):
	return meta[meta[design].isin(options)].index
	

def mutation_cohort_targets(wildcards):
    ls = []
    ls.append("analysis/fusion/merged_%s_predictions.abridged_addSample.tsv" % design),
    ls.append("analysis/fusion/%s_pyprada_fusion_table.txt" % design),
    ls.append("analysis/fusion/%s_pyprada_output.txt" % design),
    ls.append("analysis/fusion/%s_pyprada_unique_genelist.txt" % design),
    ls.append("analysis/fusion/%s_fusion_gene_table.txt" % design),
    #ls.append("analysis/fusion/%s_fusion_gene_plot.png" % design),
    #ls.append("analysis/fusion/%s_prada_homology.png" % design)
    return ls

rule mutation_cohort:
    input:
      mutation_cohort_targets

rule merge_starfusion:
    input:
      expand("analysis/fusion/{sample}/{sample}.fusion_predictions.abridged_addSample.tsv", sample=getsampleIDs(metadata))
    output:
      "analysis/fusion/merged_{design}_predictions.abridged_addSample.tsv"
    log:
      "logs/fusion/merge_{design}_fusion.log"
    message:
      "Merging fusion results"
    benchmark:
      "benchmarks/fusion/merge_{design}_fusion.benchmark"
    conda: "../envs/stat_perl_r.yml"
    params:
      meta = config['metasheet']
    shell:
      """cat {input} | sed '1 !{{/FusionName/d;}}' | awk '{{gsub(/#/,"",$1);print $0}}' > {output[0]} """

rule preprocess_prada:
    input:
      "analysis/fusion/merged_{design}_predictions.abridged_addSample.tsv"
    output:
      table = "analysis/fusion/{design}_pyprada_fusion_table.txt",
      uniquegene = "analysis/fusion/{design}_pyprada_unique_genelist.txt"
    log:
      "logs/fusion/{design}_preprocess_prada.log"
    message:
      "Preprocessing prada input file"
    benchmark:
      "benchmarks/fusion/{design}_preprocess_prada.benchmark"
    conda: "../envs/stat_perl_r.yml"
    params:
      outdir = "analysis/fusion/{design}_",
      path = "set +eu;source activate %s" % config['stat_root'],
      #pheno = config["design"],
      gtf = config['annotation_pyprada'],
      anno = 'analysis/fusion/pyprada_annotation.txt'
    shell:
      "{params.path}; Rscript src/mutation/preprocess_prada.R --fusion {input}  --outdir {params.outdir} "
      """ && cat {output.table} | sed 's/\\t/\\n/g' | sort | uniq > {output.uniquegene}   """
      """ && grep -f {output.uniquegene} {params.gtf} > {params.anno} """


rule run_prada:
    input:
      "analysis/fusion/{design}_pyprada_fusion_table.txt"
    output:
      "analysis/fusion/{design}_pyprada_output.txt"
    log:
      "logs/fusion/{design}_run_prada.log"
    message:
      "Running pyprada"
    benchmark:
      "benchmarks/fusion/{design}_pyprada.benchmark"
    conda: "../envs/stat_perl_r.yml"
    params:
      outdir = "analysis/fusion",
      config = "static/fusion/prada_config.txt",
      tmpout= "analysis/fusion/tmp/",
      prada_path = config["prada_path"],
      path = "set +eu;source activate %s" % config['prada_root']
    shell:
      "{params.path}"
      "&& {params.prada_path}/prada-homology -i {input} -o {output} -tmpdir {params.tmpout} -conf {params.config}"


rule fusion_plot:
    input:
      prada_input = "analysis/fusion/{design}_pyprada_output.txt",
      fusion = "analysis/fusion/merged_{design}_predictions.abridged_addSample.tsv",
      tpm_batch = "analysis/batchremoval/tpm.genesymbol.batchremoved.csv"
    output:
      fusion_table = "analysis/fusion/{design}_fusion_gene_table.txt",
      #fusion_plot = "analysis/fusion/design}_fusion_gene_plot.png",
      #prada_plot = "analysis/fusion/{design}_prada_homology.png"
    log:
      "logs/fusion/{design}_fusion_plot.log"
    message:
      "Running fusion plotting"
    benchmark:
      "benchmarks/fusion/{design}_fusion_plot.benchmark"     
    conda: "../envs/stat_perl_r.yml"
    params:
      outdir = "analysis/fusion/",
      meta = config["metasheet"], 
      path = "set +eu;source activate %s" % config['stat_root'],
      annotation = "static/fusion/cancerGeneList.tsv",
      phenotype = config['design']
    shell:
      "{params.path}; Rscript src/mutation/fusion_plot.R  --pradafusion {input.prada_input}  --meta {params.meta}  --expression {input.tpm_batch}  \
      --annot  {params.annotation}  --outdir  {params.outdir}  --phenotype {params.phenotype}  --input {input.fusion}"

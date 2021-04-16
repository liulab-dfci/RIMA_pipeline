#!/usr/bin/env python

#-------------------------------Immune Repertoire Individual -----------------------------#
_trust4_threads =32

def immune_repertoire_individual_targets(wildcards):
    ls = []
    for sample in config["samples"]:
        ls.append("analysis/trust4/%s/%s_cdr3.out" % (sample, sample))
        ls.append("analysis/trust4/%s/%s_cdr3.out.processed.txt" % (sample, sample))
        ls.append("analysis/trust4/%s/%s_TRUST4_BCR_heavy_cluster.Rdata" % (sample, sample))
    ls.append("files/immune_repertoire/TRUST4_BCR_light.Rdata")
    ls.append("files/immune_repertoire/TRUST4_BCR_heavy.Rdata")
    ls.append("files/immune_repertoire/TRUST4_TCR.Rdata")
    ls.append("files/immune_repertoire/TRUST4_BCR_heavy_cluster.Rdata")
    ls.append("files/immune_repertoire/TRUST4_BCR_heavy_clonality.Rdata")
    ls.append("files/immune_repertoire/TRUST4_BCR_heavy_SHMRatio.Rdata")
    ls.append("files/immune_repertoire/TRUST4_BCR_heavy_lib_reads_Infil.Rdata")
    ls.append("files/immune_repertoire/TRUST4_BCR_Ig_CS.Rdata")
    ls.append("files/immune_repertoire/TRUST4_TCR_clonality.Rdata")
    ls.append("files/immune_repertoire/TRUST4_TCR_lib_reads_Infil.Rdata")
    return ls

def cdr3_process_inputs(inputs):
    inputs_format = ','.join(str(i) for i in list(inputs)[0])
    return inputs_format

def star_process_inputs(inputs):
    inputs_format = ','.join(str(i) for i in list(inputs)[0])
    return inputs_format

def ssbcr_cluster_input(inputs):
    inputs_format = ','.join(str(i) for i in list(inputs)[0])
    return inputs_format

rule trust4_all:
    input:
      immune_repertoire_individual_targets

rule trust4_repertoire:
    input:
      "analysis/star/{sample}/{sample}.sorted.bam"
    output:
      "analysis/trust4/{sample}/{sample}_cdr3.out"
    log:
      "logs/trust4/{sample}.trust4_reper.log"
    params:
      bcrtcr_fa = config['trust4_reper_path'],
      IMGT_fa = config['trust4_IMGT_path'],
      threads = _trust4_threads,
      prefix = "analysis/trust4/{sample}/{sample}",
      trust4_path=config["trust4_path"]
    message:
      "Running TRUST4 on {wildcards.sample}"
    benchmark:
      "benchmarks/trust4/{sample}.trust4.benchmark"
    shell:
      "{params.trust4_path}/run-trust4 -f {params.bcrtcr_fa} --ref {params.IMGT_fa} -b {input} -t {params.threads} -o {params.prefix} && "
      "rm -f {params.prefix}*fq"

rule cdr3_preprocess:
    input:
      "analysis/trust4/{sample}/{sample}_cdr3.out"
    output:
      "analysis/trust4/{sample}/{sample}_cdr3.out.processed.txt"
    params:
      tmp="analysis/trust4/{sample}/{sample}_cdr3.out.processed.tmp",
      trust4_path=config["trust4_path"]
    shell:
      """perl {params.trust4_path}/trust-simplerep.pl {input} > {params.tmp}"""
      """ && sed -ig '1,1s/#count/count/g' {params.tmp} """
      """ && awk '{{print FILENAME}}' {params.tmp}  | awk  '{{print $3}}' FS='\t' | paste {params.tmp} - | awk -F '\t' 'NR==1{{$9="sample"}}1' OFS='\t'> {output}"""
      """ && rm {params.tmp} {params.tmp}g"""

rule ss_bcr_cluster:
    input:
      "analysis/trust4/{sample}/{sample}_cdr3.out.processed.txt"
    output:
      "analysis/trust4/{sample}/{sample}_TRUST4_BCR_heavy_cluster.Rdata"
    benchmark:
      "benchmarks/trust4/{sample}/{sample}_trust4_bcr_cluster.benchmark"
    params:
      meta = config["metasheet"],
      phenotype_col = config["immunerepertoire_clinical_phenotype"],
      out_dir = "analysis/trust4/{sample}/{sample}",
      wk_dir = "analysis/trust4/{sample}/",
      path= "set +eu;source activate %s" % config['stat_root'],
    conda: "../envs/stat_perl_r.yml"
    shell:
      #"cd {params.wk_dir} && "
      "{params.path}; Rscript src/immune_repertoire/trust4_bcr_cluster.R --cdr3 {input} --output {params.out_dir} --meta {params.meta} --clinic_col {params.phenotype_col}"


rule cdr3_process:
    input:
      cdr3_files=expand("analysis/trust4/{sample}/{sample}_cdr3.out.processed.txt", sample=config["samples"]),
      stat_files=expand("analysis/star/{sample}/{sample}.sorted.bam.stat.txt", sample=config["samples"]),
      ssbcr_cluster=expand("analysis/trust4/{sample}/{sample}_TRUST4_BCR_heavy_cluster.Rdata", sample=config["samples"])
    output:
      "files/immune_repertoire/TRUST4_BCR_light.Rdata",
      "files/immune_repertoire/TRUST4_BCR_heavy.Rdata",
      "files/immune_repertoire/TRUST4_TCR.Rdata",
      "files/immune_repertoire/TRUST4_BCR_heavy_cluster.Rdata",
      "files/immune_repertoire/TRUST4_BCR_heavy_clonality.Rdata",
      "files/immune_repertoire/TRUST4_BCR_heavy_SHMRatio.Rdata",
      "files/immune_repertoire/TRUST4_BCR_heavy_lib_reads_Infil.Rdata",
      "files/immune_repertoire/TRUST4_BCR_Ig_CS.Rdata",
      "files/immune_repertoire/TRUST4_TCR_clonality.Rdata",
      "files/immune_repertoire/TRUST4_TCR_lib_reads_Infil.Rdata"
    benchmark:
      "benchmarks/trust4/trust4_process.benchmark"
    log:
      "logs/trust4/trust4_process.log"
    conda: "../envs/stat_perl_r.yml"
    params:
      meta=config['metasheet'],
      outdir="files/immune_repertoire/",
      cdr3_process_input=lambda wildcards, input: cdr3_process_inputs({input.cdr3_files}),
      stat_input=lambda wildcards, input: star_process_inputs({input.stat_files}),
      ssbcr_cluster_input=lambda wildcards, input: ssbcr_cluster_input({input.ssbcr_cluster}),
      phenotype_col=config["immunerepertoire_clinical_phenotype"],
      path="set +eu;source activate %s" % config['stat_root'],
    shell:
      "{params.path}; Rscript src/immune_repertoire/trust4_process.R --cdr3 {params.cdr3_process_input} --clinic_col {params.phenotype_col} --meta {params.meta} --output {params.outdir} --stat {params.stat_input} --ssbcr_cluster {params.ssbcr_cluster_input}"

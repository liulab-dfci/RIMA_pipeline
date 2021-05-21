#!/usr/bin/env python

#-------------------------------Immune Repertoire Individual -----------------------------#
_trust4_threads =32

def immune_repertoire_individual_targets(wildcards):
    ls = []
    for sample in config["samples"]:
        ls.append("analysis/trust4/%s/%s_report.tsv" % (sample, sample))
        ls.append("analysis/trust4/%s/%s_report.processed.txt" % (sample, sample))
        
        ls.append("analysis/trust4/%s/%s_TRUST4_BCR_heavy.Rdata" % (sample, sample))
        ls.append("analysis/trust4/%s/%s_TRUST4_BCR_light.Rdata" % (sample, sample))
        ls.append("analysis/trust4/%s/%s_TRUST4_BCR_heavy_cluster.Rdata" % (sample, sample))
        ls.append("analysis/trust4/%s/%s_TRUST4_BCR_heavy_clonality.Rdata" % (sample, sample))
        ls.append("analysis/trust4/%s/%s_TRUST4_BCR_heavy_SHMRatio.Rdata" % (sample, sample))
        ls.append("analysis/trust4/%s/%s_TRUST4_BCR_heavy_lib_reads_Infil.Rdata" % (sample, sample))    
        
        ls.append("analysis/trust4/%s/%s_TRUST4_TCR.Rdata" % (sample, sample))
        ls.append("analysis/trust4/%s/%s_TRUST4_TCR_clonality.Rdata" % (sample, sample))
        ls.append("analysis/trust4/%s/%s_TRUST4_TCR_heavy_lib_reads_Infil.Rdata" % (sample, sample))
    return ls

rule trust4_all:
    input:
      immune_repertoire_individual_targets

rule trust4_repertoire:
    input:
      "analysis/star/{sample}/{sample}.sorted.bam"
    output:
      "analysis/trust4/{sample}/{sample}_report.tsv"
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
      "analysis/trust4/{sample}/{sample}_report.tsv"
    output:
      "analysis/trust4/{sample}/{sample}_report.processed.txt"
    params: 
    shell:
      """ sed '1,1s/#count/count/g' {input} > {output}"""



rule ss_bcr_process:
    input:
      report = "analysis/trust4/{sample}/{sample}_report.processed.txt",
      stat = "analysis/star/{sample}/{sample}.sorted.bam.stat.txt"
    output:
      "analysis/trust4/{sample}/{sample}_TRUST4_BCR_heavy.Rdata",
      "analysis/trust4/{sample}/{sample}_TRUST4_BCR_light.Rdata",
      "analysis/trust4/{sample}/{sample}_TRUST4_BCR_heavy_cluster.Rdata",
      "analysis/trust4/{sample}/{sample}_TRUST4_BCR_heavy_clonality.Rdata",
      "analysis/trust4/{sample}/{sample}_TRUST4_BCR_heavy_SHMRatio.Rdata",
      "analysis/trust4/{sample}/{sample}_TRUST4_BCR_heavy_lib_reads_Infil.Rdata"
    benchmark:
      "benchmarks/trust4/{sample}/{sample}_trust4_bcr_cluster.benchmark"
    params:
      out_dir = "analysis/trust4/{sample}/{sample}",
      sampleid = "{sample}",
      path= "set +eu;source activate %s" % config['stat_root'],
    conda: "../envs/stat_perl_r.yml"
    shell:
      "{params.path}; Rscript src/immune_repertoire/trust4_bcr_process.R --cdr3 {input.report} \
      --sampleid {params.sampleid} --stat {input.stat} --outdir {params.out_dir} "
      
      
rule ss_tcr_process:
    input:
      report = "analysis/trust4/{sample}/{sample}_report.processed.txt",
      stat = "analysis/star/{sample}/{sample}.sorted.bam.stat.txt"
    output:
      "analysis/trust4/{sample}/{sample}_TRUST4_TCR.Rdata",
      "analysis/trust4/{sample}/{sample}_TRUST4_TCR_clonality.Rdata",
      "analysis/trust4/{sample}/{sample}_TRUST4_TCR_heavy_lib_reads_Infil.Rdata"
    benchmark:
      "benchmarks/trust4/{sample}/{sample}_trust4_bcr_cluster.benchmark"
    params:
      out_dir = "analysis/trust4/{sample}/{sample}",
      sampleid = "{sample}",
      path= "set +eu;source activate %s" % config['stat_root'],
    conda: "../envs/stat_perl_r.yml"
    shell:
      "{params.path}; Rscript src/immune_repertoire/trust4_tcr_process.R --cdr3 {input.report} \
      --sampleid {params.sampleid} --stat {input.stat} --outdir {params.out_dir} "
      



#!/usr/bin/env python

#-------------------------------RSEM gene quantification------------------------#
# @author: Jin Wang @ref:VIPER
# @email: jwang0611@gmail.com
# @date: May, 2019
_salmon_threads = 16
import pandas as pd

metadata = pd.read_csv(config["metasheet"], index_col=0, sep=',')

def merge_sep_inputs(inputs):
    inputs_format = ' -f '.join(str(i) for i in list(inputs)[0])
    return inputs_format   

def salmon_targets(wildcards):
    ls = []
    for sample in config["samples"]:
      ls.append("analysis/salmon/%s/%s.quant.sf" % (sample, sample))
      ls.append("analysis/salmon/salmon_tpm_matrix.csv")
    return ls


rule salmon_all:
    input:
        salmon_targets

rule salmon_quantification:
    input:
        "analysis/star/{sample}/{sample}.transcriptome.bam"
    output:
        "analysis/salmon/{sample}/{sample}.quant.sf"
    log:
        "analysis/salmon/{sample}/{sample}.transcriptome.bam.log"    
    params:
        index=config["salmon_index"],
        output_path="analysis/salmon/{sample}/",
    threads: _salmon_threads
    log: "logs/salmon/{sample}.salmon_quant.log"
    #conda: "../envs/salmon_env.yml"  ##no need to use salmon in conda as salmon is installed in the common environment
    message: "salmon: from bam to sf "
    benchmark:
        "benchmarks/salmon/{sample}.salmon.benchmark"   
    shell:
        "salmon quant -t {params.index} -l A -a {input} -o {params.output_path} "
        "-p {threads} "
        " && mv {params.output_path}/quant.sf {output}"

rule salmon_matrix:
    input:
        salmon_tpm_files = expand("analysis/salmon/{sample}/{sample}.quant.sf", sample = config['samples']),
        metasheet = config["metasheet"]
    output:
        "analysis/salmon/salmon_tpm_matrix.csv"
    benchmark:
        "benchmarks/salmon/salmon_gene_matrix.benchmark"
    params:
        args = lambda wildcards, input: merge_sep_inputs({input.salmon_tpm_files})
    #conda: "../envs/salmon_env.yml" ##no need to use salmon in conda as salmon is installed in the common environment
    message: "Merge Salmon gene quantification together for all samples "
    shell:
        "perl src/raw_count_fpkm_tpm_matrix.pl --metasheet {input.metasheet} -s --header -f {params.args} 1>{output}"




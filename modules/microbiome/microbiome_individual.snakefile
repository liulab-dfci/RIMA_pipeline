#!/usr/bin/env python

#-------------------------------Microbiota Individual Classification---------------------------#
_microbiome_threads= 32

import pandas as pd

def getFastq(wildcards):
    return config["samples"][wildcards.sample]


def microbiome_individual_targets(wildcards):
    ls = []
    for sample in config["samples"]:
        ls.append("analysis/microbiome/%s/%s_classification.txt.gz" % (sample, sample))
        ls.append("analysis/microbiome/%s/%s_report.txt" % (sample, sample))
        ls.append("analysis/microbiome/%s/%s_addSample_report.txt" % (sample, sample))
    return ls

rule microbiome_individual_all:
    input:
        microbiome_individual_targets

rule centrifuge_microbiota:
    input:
        getFastq
    output:
        classfication = "analysis/microbiome/{sample}/{sample}_classification.txt.gz",
        report = "analysis/microbiome/{sample}/{sample}_report.txt",
        add_sample = "analysis/microbiome/{sample}/{sample}_addSample_report.txt"
    log:
        "logs/microbiome/{sample}.centrifuge.log"
    params:
        centrifuge_index = config["centrifuge_index"],
        path="set +eu;source activate %s" % config['centrifuge_root'],
    message:
        "Running Centrifuge on {wildcards.sample}"
    benchmark:
        "benchmarks/microbiome/{sample}.microbiome.benchmark"
    threads: _microbiome_threads     
    conda: "../envs/centrifuge_env.yml"
    shell:
        """{params.path}; centrifuge -x {params.centrifuge_index} -p {threads}  --host-taxids 9606 -1 {input[0]} -2 {input[1]} -S {output.classfication} --report-file {output.report} """
        """ && gzip {output.classfication} """
        """ && awk '{{print FILENAME}}' {output.report} | awk -F '/' '{{print $3}}' | paste - {output.report} | awk -F '\t' 'NR==1{{$1="sample"}}1' OFS='\t' > {output.add_sample} """

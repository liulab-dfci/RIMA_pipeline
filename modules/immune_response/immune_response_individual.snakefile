#!/usr/bin/env python

#-------------------------------immune response indiviudal run-------------------------#
#----------------------------MSIsensor2--------------------#
def runsMSisensorHelper(wildcards):
    r = config['runs'][wildcards.run]
    bam=["analysis/star/{sample}/{sample}.sorted.bam".format(sample=sample) for sample in r]
    index=["analysis/star/{sample}/{sample}.sorted.bam.bai".format(sample=sample) for sample in r]
    samples=bam+index
    return samples


def getSingle(wildcards):
        Single=runsMSisensorHelper(wildcards)
        if len(Single)<4:
            return Single

def msisensor_input(wildcards):
    ls=[]
    for run in config["runs"]:
        ls.append("analysis/msisensor/single/%s/%s_msisensor" % (run, run))
    return ls

def get_runs_tumor(wildcards):
    ls = []
    for run in config["runs"]:
        ls.append(config["runs"][run][0])
    return ls

def immune_response_individual_targets(wildcards):
    ls=[]
    for run in config["runs"]:
            ls.append("analysis/msisensor/single/%s/%s_msisensor" % (run, run))
            ls.append("analysis/msisensor/single/%s/%s_msisensor_dis" % (run, run))
            ls.append("analysis/msisensor/single/%s/%s_msisensor_somatic" % (run, run))
    return ls



rule immune_response_individual_all:
    input:
        immune_response_individual_targets


rule msisensor_tumor_call:
    input:
        getSingle
    output:
        "analysis/msisensor/single/{run}/{run}_msisensor",
        "analysis/msisensor/single/{run}/{run}_msisensor_dis",
        "analysis/msisensor/single/{run}/{run}_msisensor_somatic"
    message:
        "Running msisensor on {wildcards.run}"
    benchmark:
        "benchmarks/msisensor/single/{run}.msisensor.benchmark"
    log:
        "logs/msisensor/single/{run}.msisensor.log"
    params:
        prefix = "analysis/msisensor/single/{run}/{run}_msisensor",
        msisensor2_path= config["msisensor2_path"],
    shell:
        "{params.msisensor2_path}/msisensor2 msi -M {params.msisensor2_path}/models_hg38 -t {input[0]} -o {params.prefix}"

#!/usr/bin/env python

#-------------------------------MSISENSOR-------------------------
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

def msisensor_targets(wildcards):
    ls=[]
    for run in config["runs"]: 
            ls.append("analysis/msisensor/single/%s/%s_msisensor" % (run, run))
            ls.append("analysis/msisensor/single/%s/%s_msisensor_dis" % (run, run))
            ls.append("analysis/msisensor/single/%s/%s_msisensor_somatic" % (run, run))
    ls.append("files/msisensor/msi_score.txt")
    ls.append("images/msisensor/msi_score_density_plot.png")
    ls.append("images/msisensor/msi_score_comparison.png")
    return ls



rule msisensor_all:
    input:
        msisensor_targets


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


rule msisensor_plot:
    input:
        msisensor_input
    output:
        msi_score = "files/msisensor/msi_score.txt",
        msi_density = "images/msisensor/msi_score_density_plot.png",
        msi_comparison = "images/msisensor/msi_score_comparison.png"
    log:
        "analysis/variant/variant_missensor.log"
    message: 
        "Running msisensor ploting"
    benchmark:
        "analysis/variant/msisensor_plot.benchmark.txt"
    conda: "../envs/stat_perl_r.yml"
    params:
        outpath = "images/msisensor/",
        run = get_runs_tumor,
        phenotype = lambda wildcards: ','.join(str(i) for i in config["designs"]),
        meta = config["metasheet"],
        path="set +eu;source activate %s" % config['stat_root'],
    shell:  
        """cat {input} | sed '/Total_Number_of_Sites/d' | awk -F '\t' '{{print $3}}' > tmp """ 
        """&&  echo {params.run}| xargs -n1 | paste - tmp > {output.msi_score}"""
        """&& rm tmp """
        """&& {params.path};Rscript src/msi_plot.R --msiscore {output.msi_score}  --outdir {params.outpath} --phenotype {params.phenotype} --meta {params.meta}"""

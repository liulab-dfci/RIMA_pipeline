#!/usr/bin/env python

#-------------------------------STAR Fusion individual-----------------------#
metadata = pd.read_csv(config["metasheet"], index_col=0, sep=',')

def mutation_individual_targets(wildcards):
    ls = []
    for sample in config["samples"]:
        ls.append("analysis/fusion/%s/%s.fusion_predictions.tsv" % (sample, sample))
        ls.append("analysis/fusion/%s/%s.fusion_predictions.abridged.tsv" % (sample, sample))
        ls.append("analysis/fusion/%s/%s.fusion_predictions.abridged_addSample.tsv" % (sample, sample))
    return ls

def getFastq(wildcards):
    return config["samples"][wildcards.sample]

rule mutation_individual:
    input:
        mutation_individual_targets

rule STAR_fusion:
    input:
        "analysis/star/{sample}/{sample}.Chimeric.out.junction"
    output:
        res = "analysis/fusion/{sample}/{sample}.fusion_predictions.tsv",
        filter_res = "analysis/fusion/{sample}/{sample}.fusion_predictions.abridged.tsv",
        processed_res = "analysis/fusion/{sample}/{sample}.fusion_predictions.abridged_addSample.tsv"
    log:
        "analysis/fusion/{sample}/{sample}.star_fusion.log"
    shadow: "shallow"
    params:
        prefix=lambda wildcards: "analysis/fusion/{sample}".format(sample=wildcards.sample)
    message: "Running STAR fusion on {wildcards.sample}"
    benchmark:
        "analysis/fusion/{sample}/{sample}.star_fusion.benchmark.txt"
    conda: "../envs/fusion_env.yml"
    shell:
        "STAR-Fusion --chimeric_junction {input} --genome_lib_dir {config[star_fusion_index]} --output_dir {params.prefix} "
        " && mv analysis/fusion/{wildcards.sample}/star-fusion.fusion_predictions.tsv {output.res}"
        " && mv analysis/fusion/{wildcards.sample}/star-fusion.fusion_predictions.abridged.tsv {output.filter_res}"
        " && touch {output[1]}" # For some sample, final.abridged is created but not .final file; temp hack before further investigate into this
        """ && awk '{{print ARGV[1]}}' {output.filter_res} | awk -F '/' '{{print $3}}' | paste {output.filter_res} - | awk -F '\t' 'NR==1{{$16="sample"}}1' > {output.processed_res}""" ###add sample name as one column in fusion result

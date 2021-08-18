#!/usr/bin/env python

#-------------------------------PathSeq for Microbiota Classification---------------------------#
# @author: Jin Wang; @ref: VIPER
# @email: jwang0611@gmail.com

import pandas as pd

def getFastq(wildcards):
    return config["samples"][wildcards.sample]

def plot_input(wildcards):
    ls = []
    for sample in config["samples"]:
        if config["centrifuge"]:
            ls.append("analysis/centrifuge/%s/%s_addSample_report.txt" % (sample, sample))
        else:
            ls.append("analysis/pathseq/%s/%s_addSample.pathseq.txt" % (sample, sample))
    return ls

def microbe_targets(wildcards):
    ls = []
    if config["centrifuge"]:
        for sample in config["samples"]:
            ls.append("analysis/centrifuge/%s/%s_classification.txt.gz" % (sample, sample))
            ls.append("analysis/centrifuge/%s/%s_report.txt" % (sample, sample))
            ls.append("analysis/centrifuge/%s/%s_addSample_report.txt" % (sample, sample))
    else:
        for sample in config["samples"]:
            ls.append("analysis/pathseq/%s/%s_revertsam.bam" % (sample, sample))
            ls.append("analysis/pathseq/%s/%s_revertsam_addRG.bam" % (sample, sample))
            ls.append("analysis/pathseq/%s/%s_filter_metrics.txt" % (sample, sample))
            ls.append("analysis/pathseq/%s/%s_pathseq.txt" % (sample, sample))
            ls.append("analysis/pathseq/%s/%s_score_metrics.txt" % (sample, sample))
            ls.append("analysis/pathseq/%s/%s_addSample.pathseq.txt" % (sample, sample))
    ls.append("files/microbiota/merged_microbiota_abundance.txt")
    ls.append("files/microbiota/selected_microbes_taxonID.txt")
    ls.append("files/multiqc/microbiota/Microbiome_abundance_ratio.txt")
#    ls.append("files/microbiota/selected_microbes_tree.nwk")
#    ls.append("images/microbiota/selected_microbes_phylogenetic.png")
    return ls

rule microbe_all:
    input:
        microbe_targets

rule centrifuge_microbiota:
    input:
        getFastq
    output:
        classfication = "analysis/centrifuge/{sample}/{sample}_classification.txt.gz",
        report = "analysis/centrifuge/{sample}/{sample}_report.txt",
        add_sample = "analysis/centrifuge/{sample}/{sample}_addSample_report.txt"
    log:
        "logs/centrifuge/{sample}.centrifuge.log"
    params:
        # centrifuge_path = config["centrifuge_path"],
        centrifuge_index = config["centrifuge_index"],
        #mate_input = "-1 {input.fq1} -2 {input.fq2}" if len(config["samples"][metadata.index[0]]) == 2 else "-U {input.fq1} {input.fq2}"
        path="set +eu;source activate %s" % config['centrifuge_root'],
    message: 
        "Running Centrifuge on {wildcards.sample}"
    benchmark:
        "benchmarks/centrifuge/{sample}.centrifuge.benchmark"
    conda: "../envs/centrifuge_env.yml"
    shell:
        """{params.path}; centrifuge -x {params.centrifuge_index} -p 15 --host-taxids 9606 -1 {input[0]} -2 {input[1]} -S {output.classfication} --report-file {output.report} """
        """ && gzip {output.classfication} """
        """ && awk '{{print FILENAME}}' {output.report} | awk -F '/' '{{print $3}}' | paste - {output.report} | awk -F '\t' 'NR==1{{$1="sample"}}1' OFS='\t' > {output.add_sample} """

rule gatk_ubam:
    input:
        "analysis/star/{sample}/{sample}.sorted.bam"
    output:
        "analysis/pathseq/{sample}/{sample}_revertsam.bam"
    message: 
        "Running Picard RevertSam on {wildcards.sample}"
    benchmark:
        "benchmarks/pathseq/{sample}.RevertSam.benchmark"
    log:
        "logs/pathseq/{sample}.RevertSam.log"
    params:
        path="set +eu;source activate %s" % config['gatk4_root'],
    conda: "../envs/gatk4_env.yml"
    shell:
        "{params.path}; gatk RevertSam -I {input} -O {output}"

rule picard_addReadGroup:
    input:
        "analysis/pathseq/{sample}/{sample}_revertsam.bam"
    output:
        "analysis/pathseq/{sample}/{sample}_revertsam_addRG.bam"
    message: 
        "Running Picard AddReadsGroup on {wildcards.sample}"
    benchmark:
        "benchmarks/pathseq/{sample}.addReadsGroup.benchmark"
    log:
        "logs/pathseq/{sample}.addReadsGroup.log"
    params:
        rgID = "H0164.2",
        rgLB = "Solexa-272222",
        rgPL = "illumina",
        rgPU = "H0164ALXX140820.2",
        rgSM = "{wildcards.sample}",
        path="set +eu;source activate %s" % config['gatk4_root'],
    conda: "../envs/gatk4_env.yml"
    shell:
        "{params.path}; picard AddOrReplaceReadGroups I= {input} O= {output} SORT_ORDER=coordinate "   
        "RGID={params.rgID} RGLB={params.rgLB} RGPL={params.rgPL} RGPU={params.rgPU} RGSM={params.rgSM} CREATE_INDEX=True"
        " && rm {input}"

rule pathseq_microbiota_call:
    input:
        bam = "analysis/pathseq/{sample}/{sample}_revertsam_addRG.bam"
    output:
        filter_metrics ="analysis/pathseq/{sample}/{sample}_filter_metrics.txt",
        score_output = "analysis/pathseq/{sample}/{sample}_pathseq.txt",
        score_metrics = "analysis/pathseq/{sample}/{sample}_score_metrics.txt",
        add_sample="analysis/pathseq/{sample}/{sample}_addSample_pathseq.txt"
    message: 
        "Running PathSeq on {wildcards.sample}"
    params:
        min_score_ident = "0.85",
        min_clipped_read_len = "70",
        microbe_fa = config["microbe_fa_path"],
        microbe_fa_img = config["microbe_img_path"],
        taxonomy = config["taxonomy_path"],
        host_kmer = config["host_kmer_path"],
        host_fa_img = config["host_img_path"],
        path="set +eu;source activate %s" % config['gatk4_root'],
    benchmark:
        "benchmarks/pathseq/{sample}.pathseq.benchmark"
    log:
        "logs/pathseq/{sample}.pathseq.log"
    conda: "../envs/gatk4_env.yml"
    shell:
        "{params.path}; gatk PathSeqPipelineSpark "
        "--input {input.bam} "
        "--kmer-file {params.host_kmer} "
        "--filter-bwa-image {params.host_fa_img} "
        "--microbe-fasta {params.microbe_fa} "
        "--microbe-bwa-image {params.microbe_fa_img} "
        "--taxonomy-file {params.taxonomy} "
        "--scores-output {output.score_output} "
        "--min-score-identity {params.min_score_ident} "
        "--min-clipped-read-length {params.min_clipped_read_len} "
        "--filter-metrics {output.filter_metrics} "
        "--score-metrics {output.score_metrics} "
        "--create-output-bam-index false"
        """ && awk '{{print FILENAME}}' {output.score_output} | awk -F '/' '{{print $3}}' | paste - {output.score_output} | awk -F '\t' 'NR==1{{$1="sample"}}1' OFS='\t' > {output.add_sample}"""

rule micorbiota_plot:
    input:
        plot_input
    output:
        merged_file = "files/microbiota/merged_microbiota_abundance.txt",
        selected_ID = "files/microbiota/selected_microbes_taxonID.txt",
        multiqc = "files/multiqc/microbiota/Microbiome_abundance_ratio.txt"
    params:
        outdir = "images/microbiota/",
        phenotype = config['designs'],
        mul_folder = "files/multiqc/microbiota/",
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
        " && {params.path}; Rscript src/mic_plot_final.R --input {output.merged_file} --outdir {params.outdir} --clinic_col {params.phenotype} --meta {params.meta} --multiqc {params.mul_folder}"
        " && mv {params.outdir}selected_microbes_taxonID.txt {output.selected_ID}"

        
#rule taxon_to_tree:
#    input:
#        "files/microbiota/selected_microbes_taxonID.txt"
#    output:
#        "files/microbiota/selected_microbes_tree.nwk"     
#    message: 
#        "Running Microbiota taxon to tree"
#    benchmark:
#        "benchmarks/centrifuge/centrifuge_taxon_tree.benchmark"
#    params:
#        path="set +eu;source activate %s" % config['centrifuge_root'],
#    conda: "../envs/centrifuge_env.yml"
#    shell:
#        "cut -f1 {input} | {params.path}; ete3 ncbiquery --tree > {output} "
#
#rule tree_plot:
#    input:
#        "files/microbiota/selected_microbes_tree.nwk"
#    output:
#        "images/microbiota/selected_microbes_phylogenetic.png"
#    params:
#        outdir = "images/microbiota/",
#        meta = config["metasheet"],
#        path="set +eu;source activate %s" % config['stat_root'],
#    message: 
#        "Running Microbiota tree plot"
#    benchmark:
#        "benchmarks/centrifuge/centrifuge_tree_plot.benchmark"
#    conda: "../envs/stat_perl_r.yml"
#    shell:
#        "{params.path}; Rscript src/tree_plot.R --input {input} --outdir {params.outdir} "
#

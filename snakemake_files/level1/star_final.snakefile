#!/usr/bin/env python

#-------------------------------STAR Alignment------------------------#
# @author: Jin Wang @ref:VIPER
# @email: jwang0611@gmail.com

import pandas as pd

metadata = pd.read_csv(config["metasheet"], index_col=0, sep=',')

###input data(compressed or not)
gz_command = "--readFilesCommand zcat" if config["samples"][metadata.index[0]][0][-3:] == '.gz' else ""

###extract gene reads count column from star output
if config["library_type"] == "fr-firststrand":
    count_col = 3
elif config["library_type"] == "fr-secondstrand":
    count_col = 2
else:
    count_col = 1

def merge_sep_inputs(inputs):
    inputs_format = ' -f '.join(str(i) for i in list(inputs)[0])
    return inputs_format

def align_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for sample in config["samples"]:
        ls.append("analysis/star/%s/%s.unsorted.bam" % (sample,sample))
        ls.append("analysis/star/%s/%s.sorted.bam" % (sample, sample))
        ls.append("analysis/star/%s/%s.transcriptome.bam" % (sample, sample))
        ls.append("analysis/star/%s/%s.Chimeric.out.junction" % (sample, sample))
        ls.append("analysis/star/%s/%s.Log.final.out" % (sample, sample))
        ls.append("analysis/star/%s/%s.counts.tab" % (sample, sample))
        ls.append("analysis/star/%s/%s.sorted.bam.stat.txt" % (sample, sample))
        ls.append("analysis/star/%s/%s.sorted.bam.bai" % (sample, sample))
    return ls

def star_report_tagets(wildcards):
    ls = []
    ls.append("files/star/STAR_Align_Report.csv" )
    ls.append("images/star/STAR_Align_Report.png" )
    ls.append("analysis/star/STAR_Gene_Counts.csv")
    return ls


def getFastq(wildcards):
    return config["samples"][wildcards.sample]


rule align_all:
    input:
        align_targets

rule report_all:
    input:
        star_report_tagets


rule star_align:
    """star alignment for RNA-Seq raw data"""
    input:
        getFastq
    output:
        unsortedBAM = "analysis/star/{sample}/{sample}.unsorted.bam",
        sortedBAM = "analysis/star/{sample}/{sample}.sorted.bam",
        transcriptomeBAM = "analysis/star/{sample}/{sample}.transcriptome.bam",
        junction_file = "analysis/star/{sample}/{sample}.Chimeric.out.junction",
        counts = "analysis/star/{sample}/{sample}.counts.tab",
        log_file = "analysis/star/{sample}/{sample}.Log.final.out"
    params:
        gz_support = gz_command,
        prefix = lambda wildcards: "analysis/star/{sample}/{sample}".format(sample=wildcards.sample)
    message:
        "Running STAR Alignment on {wildcards.sample}"
    log:
        "logs/star/{sample}.star_align.log"
    benchmark:
        "benchmarks/star/{sample}.star_align.benchmark"
    #conda: "../envs/fusion_env.yml"
    shell:
        "STAR --runThreadN {config[threads]} --genomeDir {config[star_index]} "
        "--outReadsUnmapped None "
        "--chimSegmentMin 12 "
        "--chimJunctionOverhangMin 12 "
        "--chimOutJunctionFormat 1 "
        "--alignSJDBoverhangMin 10 "
	"--alignMatesGapMax 1000000 "
        "--alignIntronMax 1000000 "
        "--alignSJstitchMismatchNmax 5 -1 5 5 "
        "--outSAMstrandField intronMotif "
        "--outSAMunmapped Within "
        "--outSAMtype BAM Unsorted "
        "--readFilesIn {input} "
        "--chimMultimapScoreRange 10 "
        "--chimMultimapNmax 10 "
        "--chimNonchimScoreDropMin 10 "
        "--peOverlapNbasesMin 12 "
        "--peOverlapMMp 0.1 "
        "--genomeLoad NoSharedMemory "
        "--outSAMheaderHD @HD VN:1.4 "
        "--twopassMode Basic "
        "{params.gz_support} "
        "--outFileNamePrefix {params.prefix} "
        "--quantMode TranscriptomeSAM GeneCounts"
        " && mv {params.prefix}Aligned.out.bam {output.unsortedBAM}"
        " && samtools sort -T {params.prefix}TMP -o {output.sortedBAM} -@ 8  {output.unsortedBAM} "
        " && mv {params.prefix}Aligned.toTranscriptome.out.bam {output.transcriptomeBAM}"
        " && mv {params.prefix}ReadsPerGene.out.tab {output.counts}"
        " && mv {params.prefix}Chimeric.out.junction {output.junction_file}"
        " && mv {params.prefix}Log.final.out {output.log_file}"


rule index_bam:
    """INDEX the {sample}.sorted.bam file"""
    input:
        "analysis/star/{sample}/{sample}.sorted.bam"
    output:
        "analysis/star/{sample}/{sample}.sorted.bam.bai"
    message:
        "Indexing {wildcards.sample}.sorted.bam"
    benchmark:
        "benchmarks/star/{sample}.index.benchmark"
    shell:
        "samtools index {input} > {output}"


rule align_bam_stat:
    input:
        bam = "analysis/star/{sample}/{sample}.sorted.bam",
        bai = "analysis/star/{sample}/{sample}.sorted.bam.bai"
    output:
        "analysis/star/{sample}/{sample}.sorted.bam.stat.txt"
    message:
        "Generating Aligned BAM stats"
    benchmark:
        "benchmarks/star/{sample}_bam_stat.benchmark"
    shell:
        "samtools stats {input.bam}| grep ^SN | cut -f 2- > {output}"



###STAR reports takes all the samples and generates the figures.Used only if want to show all the samples together
rule generate_STAR_report:
    input:
        star_log_files=expand( "analysis/star/{sample}/{sample}.Log.final.out", sample=config["samples"] ),
        star_gene_count_files=expand( "analysis/star/{sample}/{sample}.counts.tab", sample=config["samples"] )
    output:
        csv="files/star/STAR_Align_Report.csv",
        png="images/star/STAR_Align_Report.png",
        gene_counts="analysis/star/STAR_Gene_Counts.csv"
    message:
        "Generating STAR report"
    benchmark:
        "benchmarks/star/generate_STAR_report.benchmark"
    conda: "../envs/stat_perl_r.yml"
    params:
        extract_col = count_col,
        log_files = lambda wildcards, input: merge_sep_inputs({input.star_log_files}),
        count_files = lambda wildcards, input: merge_sep_inputs({input.star_gene_count_files}),
        path="set +eu;source activate %s" % config['stat_root'],
    shell:
        """{params.path}; perl src/STAR_reports.pl -f {params.log_files} 1>{output.csv} && """
        """{params.path}; Rscript src/map_stats.R -i {output.csv} -o {output.png} && """
        """{params.path}; perl src/raw_count_fpkm_tpm_matrix.pl --column {params.extract_col} -f {params.count_files} 1>{output.gene_counts}"""






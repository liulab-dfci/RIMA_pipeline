#!/usr/bin/env python

#-------------------------------Preprocessing individual module------------------------#
#---------------------------Preprocess individual targets----------------------------------#
_preprocess_threads = 64

import pandas as pd
metadata = pd.read_csv(config["metasheet"], index_col=0, sep=',')

###input data(compressed or not)
gz_command = "--readFilesCommand zcat" if config["samples"][metadata.index[0]][0][-3:] == '.gz' else ""
###rseqc ref (housekeeping or not)
rseqc_ref = config["housekeeping_bed_path"] if config["rseqc_ref"] == "house_keeping" else config["bed_path"]


def merge_sep_inputs(inputs):
    inputs_format = ' -f '.join(str(i) for i in list(inputs)[0])
    return inputs_format


def DownsamplingOrNot(bam_stat):
    sequences = open(bam_stat,"r").readlines()[2]
    seq_num = int(sequences.split(":")[1].replace("\t","").replace("\n",""))/1000000
    return seq_num



def preprocess_individual_targets(wildcards):
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
        ls.append("analysis/salmon/%s/%s.quant.sf" % (sample, sample))
        ls.append("analysis/rseqc/%s/%s.stat_tmp.txt" % (sample, sample))
        ls.append("analysis/rseqc/%s/%s_downsampling_housekeeping.bam" % (sample, sample))
        ls.append("analysis/rseqc/%s/%s_downsampling_housekeeping.bam.bai" % (sample, sample))
        ls.append("analysis/rseqc/read_distrib/%s/%s.txt" % (sample, sample))
        ls.append("files/rseqc/gene_body_cvg/%s/%s.geneBodyCoverage.r" % (sample, sample))
        ls.append("files/rseqc/junction_saturation/%s/%s.junctionSaturation_plot.pdf" % (sample, sample))
        ls.append("analysis/rseqc/insert_size/%s/metrics/%s.insert_size_metrics.txt" % (sample, sample))
        ls.append("analysis/rseqc/tin_score/%s/%s.summary.txt" % (sample, sample))
        ls.append("analysis/rseqc/tin_score/%s/%s.tin_score.txt" % (sample, sample))
    return ls

def getFastq(wildcards):
    return config["samples"][wildcards.sample]

def DownsamplingOrNot(bam_stat):
    sequences = open(bam_stat,"r").readlines()[2]
    seq_num = int(sequences.split(":")[1].replace("\t","").replace("\n",""))/1000000
    return seq_num

def getHousekeepingBam(hk):
    if hk == "house_keeping":
      return "analysis/star/{sample}/{sample}_downsampling_housekeeping.bam"
    else:
      return "analysis/star/{sample}/{sample}_downsampling.bam"

def getHousekeepingBai(hk):
    if hk == "house_keeping":
       return "analysis/star/{sample}/{sample}_downsampling_housekeeping.bam.bai"
    else:
      return "analysis/star/{sample}/{sample}_downsampling.bam.bai"

#---------------------STAR alignment rules-------------------------#
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
    threads: _preprocess_threads
    message:
        "Running STAR Alignment on {wildcards.sample}"
    log:
        "logs/star/{sample}.star_align.log"
    benchmark:
        "benchmarks/star/{sample}.star_align.benchmark"
    conda:
        "../envs/star_env.yml"
    shell:
        "STAR --runThreadN {threads} --genomeDir {config[star_index]} "
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
    conda:
        "../envs/star_env.yml"
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
    conda:
        "../envs/star_env.yml"
    shell:
        "samtools stats {input.bam}| grep ^SN | cut -f 2- > {output}"


#---------------------RSeQC rules--------------------------#

rule Size_downsampling:
    input:
        "analysis/star/{sample}/{sample}.sorted.bam.stat.txt"
    output:
        "analysis/rseqc/{sample}/{sample}.stat_tmp.txt"
    shell:
        """chmod +x src/preprocess/ds_check_size.sh && """
        """src/preprocess/ds_check_size.sh {input} {output} """

rule bam_downsampling:
    input:
        bam = "analysis/star/{sample}/{sample}.sorted.bam",
        stat = "analysis/star/{sample}/{sample}.sorted.bam.stat.txt",
        stat_tmp = "analysis/rseqc/{sample}/{sample}.stat_tmp.txt"
    output:
        Downsampling_bam = "analysis/star/{sample}/{sample}_downsampling.bam",
        Downsampling_bai = "analysis/star/{sample}/{sample}_downsampling.bam.bai"
    message:
        "Running Downsampling on {wildcards.sample}"
    log:
        "logs/star/{sample}.downsampling.log"
    benchmark:
        "benchmarks/star/{sample}.downsampling.benchmark"
    params:
        prefix = "analysis/star/{sample}",
        path="set +eu;source activate %s" % config['rseqc_root'],
    conda: "../envs/rseqc_env.yml"
    shell:
        """size=$(python -c "print(open('{input.stat_tmp}','r').readlines()[0].replace('\\n',''))") && """
        """chmod +x src/preprocess/downsampling.sh && """
        """src/preprocess/downsampling.sh {input.bam} $size && """
        """mv {input.bam}*downsampling.bam {output.Downsampling_bam} && """
        """samtools index {output.Downsampling_bam} > {output.Downsampling_bai}"""

rule Downsampling_HouseKeeping:
    input:
        bam = "analysis/star/{sample}/{sample}_downsampling.bam",
        stat_tmp = "analysis/rseqc/{sample}/{sample}.stat_tmp.txt"
    output:
        downsampling_hp_bam = "analysis/rseqc/{sample}/{sample}_downsampling_housekeeping.bam",
        Downsampling_hp_bai = "analysis/rseqc/{sample}/{sample}_downsampling_housekeeping.bam.bai"
    message:
        "Running Downsampling on house keeping genes"
    params:
        housekeeping_bed = config["housekeeping_bed_path"]
    log:
        "logs/rseqc/{sample}.downsampling_housekeeping.log"
    benchmark:
        "benchmarks/rseqc/{sample}.downsampling_housekeeping.benchmark"
    shell:
        "bedtools intersect -a {input.bam} -b {params.housekeeping_bed} > {output.downsampling_hp_bam} && "
        "samtools index {output.downsampling_hp_bam} > {output.Downsampling_hp_bai}"

     
rule tin_score:
    input:
        bam = getHousekeepingBam,
        bai = getHousekeepingBai
    output:
        summary = "analysis/rseqc/tin_score/{sample}/{sample}.summary.txt",
        score = "analysis/rseqc/tin_score/{sample}/{sample}.tin_score.txt"
    message:
        "Running RseQC TIN score on {wildcards.sample}"
    log:
        "logs/rseqc/tin_score/{sample}.tin_score.log"
    benchmark:
        "benchmarks/rseqc/tin_score/{sample}.tin_score.benchmark"
    params:
        bed_ref = rseqc_ref,
        prefix = "./{sample}",
        min_coverage = "10" ,
        sample_size = "100" ,
        path="set +eu;source activate %s" % config['rseqc_root'],
    conda: "../envs/rseqc_env.yml"
    shell:
        "{params.path}; tin.py"
        " --input={input.bam}"
        " --refgene={params.bed_ref}"
        " --minCov={params.min_coverage}"
        " --sample-size={params.sample_size} "
        "&& mv {params.prefix}*summary.txt {output.summary} "
        "&& mv {params.prefix}*tin.xls {output.score}"


rule read_distrib_qc:
    input:
        bam = getHousekeepingBam,
        bai = getHousekeepingBai
    output:
        "analysis/rseqc/read_distrib/{sample}/{sample}.txt"
    message:
        "Running RseQC read distribution on {wildcards.sample}"
    log:
        "logs/rseqc/read_distrib/{sample}.read_distrib_qc_matrix.log"
    benchmark:
        "benchmarks/rseqc/read_distrib/{sample}.read_distrib_qc_matrix.benchmark"
    params:
        bed_ref = rseqc_ref,
        path="set +eu;source activate %s" % config['rseqc_root'],
    conda: "../envs/rseqc_env.yml"
    shell:
        "{params.path}; read_distribution.py"
        " --input-file={input.bam}"
        " --refgene={params.bed_ref} 1>{output}"


rule gene_body_cvg_qc:
    input:
        bam = getHousekeepingBam,
        bai = getHousekeepingBai
    output:
        "files/rseqc/gene_body_cvg/{sample}/{sample}.geneBodyCoverage.r"
    threads: _preprocess_threads
    message:
        "Creating gene body coverage curves"
    log:
        "logs/rseqc/gene_body_cvg/{sample}.gene_body_cvg_qc.log"
    benchmark:
        "benchmarks/rseqc/gene_body_cvg/{sample}.gene_body_cvg_qc.benchmark"
    params:
        bed_ref = rseqc_ref,
        prefix = "files/rseqc/gene_body_cvg/{sample}/{sample}",
        path="set +eu;source activate %s" % config['rseqc_root'],
    conda: "../envs/rseqc_env.yml"
    shell:
        "{params.path}; geneBody_coverage.py -i {input.bam} -r {params.bed_ref}"
        " -f png -o {params.prefix}"

rule junction_saturation:
    input:
        bam = getHousekeepingBam,
        bai = getHousekeepingBai
    output:
        "files/rseqc/junction_saturation/{sample}/{sample}.junctionSaturation_plot.pdf"
    message:
        "Determining junction saturation for {wildcards.sample}"
    benchmark:
        "benchmarks/rseqc/junction_saturation/{sample}.junction_saturation.benchmark"
    log:
        "logs/rseqc/junction_saturation/{sample}.junction_saturation.log"
    params:
        prefix = "files/rseqc/junction_saturation/{sample}/{sample}",
        path="set +eu;source activate %s" % config['rseqc_root'],
    conda: "../envs/rseqc_env.yml"
    shell:
        "{params.path}; junction_saturation.py -i {input.bam} -r {config[bed_path]} -o {params.prefix}"

rule collect_insert_size:
    input:
        bam = getHousekeepingBam,
        bai = getHousekeepingBai
    output:
        text="analysis/rseqc/insert_size/{sample}/metrics/{sample}.insert_size_metrics.txt",
        pdf="analysis/rseqc/insert_size/{sample}/hist/{sample}.insert_size_histogram.pdf"
    message:
        "Collecting insert size for {wildcards.sample}"
    benchmark:
        "benchmarks/rseqc/insert_size/{sample}.collect_insert_size.benchmark"
    log:
        "logs/rseqc/insert_size/{sample}.collect_insert_size.log"
    params:
        path="set +eu;source activate %s" % config['rseqc_root'],
    conda: "../envs/rseqc_env.yml"
    shell:
        "{params.path}; picard CollectInsertSizeMetrics O={output.text} I={input.bam}  R={config[fasta_path]} M=0.5 H={output.pdf}"

#----------------------Salmon rules----------------------#
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
    threads: _preprocess_threads
    log: "logs/salmon/{sample}.salmon_quant.log"
    message: "salmon: from bam to sf "
    benchmark:
        "benchmarks/salmon/{sample}.salmon.benchmark"
    shell:
        "salmon quant -t {params.index} -l A -a {input} -o {params.output_path} "
        "-p {threads} "
        " && mv {params.output_path}/quant.sf {output}"

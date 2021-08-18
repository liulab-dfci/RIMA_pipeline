#!/usr/bin/env python

#-------------------------------
# @author: Jin Wang; @ref: VIPER
# @email: jwang0611@gmail.com

###select rseqc reference file(downsampling or not)
rseqc_ref = config["housekeeping_bed_path"] if config["rseqc_ref"] == "house_keeping" else config["bed_path"]


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

def merge_sep_inputs(inputs):
    inputs_format = ' -f '.join(str(i) for i in list(inputs)[0])
    return inputs_format   


def rseqc_targets(wildcards):
    ls = []
    for sample in config["samples"]:
        ls.append("analysis/rseqc/%s/%s.stat_tmp.txt" % (sample, sample))
        ls.append("analysis/rseqc/%s/%s_downsampling_housekeeping.bam" % (sample, sample))
        ls.append("analysis/rseqc/%s/%s_downsampling_housekeeping.bam.bai" % (sample, sample))       
        ls.append("analysis/rseqc/read_distrib/%s/%s.txt" % (sample, sample))
        ls.append("images/rseqc/gene_body_cvg/%s/%s.geneBodyCoverage.curves.png" % (sample, sample))
        ls.append("images/rseqc/gene_body_cvg/%s/%s.geneBodyCoverage.r" % (sample, sample))
        ls.append("images/rseqc/junction_saturation/%s/%s.junctionSaturation_plot.pdf" % (sample, sample))
        ls.append("images/rseqc/insert_size/%s/%s.histogram.pdf" % (sample, sample))
        ls.append("analysis/rseqc/tin_score/%s/%s.summary.txt" % (sample, sample)) 
        ls.append("analysis/rseqc/tin_score/%s/%s.tin_score.txt" % (sample, sample))  
        ls.append("images/rseqc/read_quality/%s/%s.qual.boxplot.pdf" %(sample, sample))
        ls.append("images/rseqc/read_quality/%s/%s.qual.heatmap.pdf" %(sample, sample))       
    ls.append("analysis/rseqc/gene_body_cvg/geneBodyCoverage.r")
    #ls.append("images/rseqc/gene_body_cvg/geneBodyCoverage.heatMap.png")
    ls.append("images/rseqc/gene_body_cvg/geneBodyCoverage.curves.png")
    ls.append("files/rseqc/tin_score/tin_score_summary.txt")
    ls.append("images/rseqc/tin_score/medTIN_score_plot.png")
    ls.append("analysis/rseqc/read_distrib/read_distrib.matrix.tab")
    ls.append("images/rseqc/read_distrib/read_distrib.pdf")

        
    return ls

rule rseqc_all:
    input:
        rseqc_targets

rule Size_downsampling:
    input:
        "analysis/star/{sample}/{sample}.sorted.bam.stat.txt"
    output:
        "analysis/rseqc/{sample}/{sample}.stat_tmp.txt"   
    shell:
        """chmod +x src/ds_check_size.sh && """
        """src/ds_check_size.sh {input} {output} """

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
    shell:
        """size=$(python -c "print(open('{input.stat_tmp}','r').readlines()[0].replace('\\n',''))") && """
        """chmod +x src/downsampling.sh && """
        """src/downsampling.sh {input.bam} $size && """
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
        min_coverage = "10" , ##why 10 ? Find the reason for this for Yang and Aashna
        sample_size = "100" , ##why 100 ? Find the reason for this for Yang and Aashna
        path="set +eu;source activate %s" % config['rseqc_root'],
    shell:
        "{params.path}; tin.py"
        " --input={input.bam}"
        " --refgene={params.bed_ref}"
        " --minCov={params.min_coverage}"
        " --sample-size={params.sample_size} "
        "&& mv {params.prefix}*summary.txt {output.summary} "
        "&& mv {params.prefix}*tin.xls {output.score}"

rule tin_summary:
    input:
        expand("analysis/rseqc/tin_score/{sample}/{sample}.summary.txt", sample=config['samples'])
    output:
        score = "files/rseqc/tin_score/tin_score_summary.txt",
        tin_plot = "images/rseqc/tin_score/medTIN_score_plot.png"
    message: 
        "plotting TIN score summary"
    log:
        "logs/rseqc/tin_score/tin_score_summary.log"
    benchmark:
        "benchmarks/rseqc/tin_score/tin_score_summary.benchmark"
    params: 
        outpath = "images/rseqc/tin_score/",
        path="set +eu;source activate %s" % config['stat_root'],
    conda: "../envs/stat_perl_r.yml"
    shell:
        """cat {input} | sed '1 !{{/Bam_file/d;}}' >{output.score} && """
        """{params.path}; Rscript src/tin_score_plot.R --tin_score {output.score} --outdir {params.outpath}"""

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
    shell:
        "{params.path}; read_distribution.py"
        " --input-file={input.bam}"
        " --refgene={params.bed_ref} 1>{output}"

rule read_quality:
    input:
        bam = getHousekeepingBam,
        bai = getHousekeepingBai
    output:
        "images/rseqc/read_quality/{sample}/{sample}.qual.boxplot.pdf",
        "images/rseqc/read_quality/{sample}/{sample}.qual.heatmap.pdf"        
    # priority: 4
    message: 
        "Running RseQC read quality on {wildcards.sample}"
    log:
        "logs/rseqc/read_quality/{sample}.read_quality.log"
    benchmark:
        "benchmarks/rseqc/read_quality/{sample}.read_quality.benchmark"
    params: 
        output = "images/rseqc/read_quality/{sample}/{sample}",
        path="set +eu;source activate %s" % config['rseqc_root'],
    shell:
        "{params.path}; read_quality.py -i {input.bam} -o {params.output}"


rule read_distrib_qc_matrix:
    input:
        read_distrib_files=expand( "analysis/rseqc/read_distrib/{sample}/{sample}.txt", sample = config['samples'])
    output:
        matrix="analysis/rseqc/read_distrib/read_distrib.matrix.tab",
        png="images/rseqc/read_distrib/read_distrib.pdf"
    message: 
        "Creating RseQC read distribution matrix"
    log:
        "logs/rseqc/read_distrib/read_distrib_qc_matrix.log"
    benchmark:
        "benchmarks/rseqc/read_distrib/read_distrib_qc_matrix.benchmark"
    conda: "../envs/stat_perl_r.yml"
    params:
        file_list_with_flag = lambda wildcards, input: merge_sep_inputs({input.read_distrib_files}),
        path="set +eu;source activate %s" % config['stat_root'],
    shell:
        "perl src/read_distrib_matrix.pl -f {params.file_list_with_flag} 1>{output.matrix} && "
        "{params.path}; Rscript src/read_distrib.R {output.matrix} {output.png}"
    
rule gene_body_cvg_qc:
    input:
        bam = getHousekeepingBam,
        bai = getHousekeepingBai
    output:
        "images/rseqc/gene_body_cvg/{sample}/{sample}.geneBodyCoverage.curves.png",
        "images/rseqc/gene_body_cvg/{sample}/{sample}.geneBodyCoverage.r"
    threads: 1
    message: 
        "Creating gene body coverage curves"
    log:
        "logs/rseqc/gene_body_cvg/{sample}.gene_body_cvg_qc.log"
    benchmark:
        "benchmarks/rseqc/gene_body_cvg/{sample}.gene_body_cvg_qc.benchmark"
    params: 
        bed_ref = rseqc_ref,
        prefix = "images/rseqc/gene_body_cvg/{sample}/{sample}",
        path="set +eu;source activate %s" % config['rseqc_root'],
    shell:
        "{params.path}; geneBody_coverage.py -i {input.bam} -r {params.bed_ref}"
        " -f png -o {params.prefix}"


rule plot_gene_body_cvg:
    input:
        samples_list=expand("images/rseqc/gene_body_cvg/{sample}/{sample}.geneBodyCoverage.r", sample=config["samples"] )
    output:
        rscript="analysis/rseqc/gene_body_cvg/geneBodyCoverage.r",
     #   png="images/rseqc/gene_body_cvg/geneBodyCoverage.heatMap.png",
        png_curves="images/rseqc/gene_body_cvg/geneBodyCoverage.curves.png"
    message: "Plotting gene body coverage"
    benchmark:
        "benchmarks/rseqc/gene_body_cvg/plot_gene_body_cvg.benchmark"
    params:
        path_1="set +eu;source activate %s" % config['stat_root'],
        path_2="set +eu;source activate %s" % config['rseqc_root'],
        
    conda: "../envs/stat_perl_r.yml"
    shell:
        "{params.path_2}; perl src/plot_gene_body_cvg.pl --rfile {output.rscript} --curves_png {output.png_curves}"
        " {input.samples_list} && {params.path_1}; Rscript {output.rscript}"

rule junction_saturation:
    input:
        bam = getHousekeepingBam,
        bai = getHousekeepingBai
    output:
        "images/rseqc/junction_saturation/{sample}/{sample}.junctionSaturation_plot.pdf"
    message: 
        "Determining junction saturation for {wildcards.sample}"
    benchmark:
        "benchmarks/rseqc/junction_saturation/{sample}.junction_saturation.benchmark"
    log:
        "logs/rseqc/junction_saturation/{sample}.junction_saturation.log"    
    params: 
        prefix = "images/rseqc/junction_saturation/{sample}/{sample}",
        path="set +eu;source activate %s" % config['rseqc_root'],
    shell:
        "{params.path}; junction_saturation.py -i {input.bam} -r {config[bed_path]} -o {params.prefix}"


rule collect_insert_size:
    input:
        bam = getHousekeepingBam,
        bai = getHousekeepingBai
    output:
        "images/rseqc/insert_size/{sample}/{sample}.histogram.pdf"
    message: 
        "Collecting insert size for {wildcards.sample}"
    benchmark:
        "benchmarks/rseqc/insert_size/{sample}.collect_insert_size.benchmark"
    log:
        "logs/rseqc/insert_size/{sample}.collect_insert_size.log"
    params:
        prefix = "images/rseqc/insert_size/{sample}/{sample}",
        path="set +eu;source activate %s" % config['rseqc_root'],
    shell:
        "{params.path}; picard CollectInsertSizeMetrics H={output} I={input.bam} O={params.prefix} R={config[fasta_path]}"


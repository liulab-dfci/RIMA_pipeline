configfile: "config.yaml"

def bam2fastq_targets(wildcards):
    ls = []
    for sample in config['samples']:
        ls.append("data/fastqs/%s_1.fq" % sample)
    return ls

def getBam(wildcards):
    sample = config['samples'][wildcards.sample]
    #print(sample)
    return sample

rule all:
    input:
        bam2fastq_targets

rule bam2fastq:
    input:
        getBam
    output:
        fq1="data/fastqs/{sample}_1.fq",
        fq2="data/fastqs/{sample}_2.fq",
    threads: 4
    benchmark: "benchmarks/bam2fastq_{sample}.txt"
    shell:
        "samtools view -H {input} | grep \"^@RG\" > {wildcards.sample}.header && samtools collate -@ 32 -Ouf {input} | samtools fastq -@ 32 -1 {output.fq1} -2 {output.fq2} -t -s /dev/null -0 /dev/null -"

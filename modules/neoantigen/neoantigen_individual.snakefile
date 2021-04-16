#!/usr/bin/env python

#-------------------------------Neoantigen individual module----------------#

_neoantigen_individual_threads=16

def neoantigen_individual_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for sample in config["samples"]:
        ls.append("analysis/neoantigen/%s/%s.alignment.p" % (sample,sample))
        ls.append("analysis/neoantigen/%s/%s.genes.json" % (sample,sample))
        ls.append("analysis/neoantigen/%s/%s.genotype.json" % (sample,sample))
        ls.append("analysis/neoantigen/%s/%s.genotype.log" % (sample,sample))
        ls.append("analysis/neoantigen/%s/%s.sorted.extracted.1.fq.gz" % (sample,sample))
        ls.append("analysis/neoantigen/%s/%s.sorted.extracted.2.fq.gz" % (sample,sample))
        ls.append("analysis/neoantigen/%s/%s.1.fq.gz" % (sample,sample))
        ls.append("analysis/neoantigen/%s/%s.2.fq.gz" % (sample,sample))
        ls.append("analysis/neoantigen/merge/%s.genotype.json" % (sample))
    return ls

rule neoantigen_individual_all:
     input:
         neoantigen_individual_targets

###------------------arcasHLA individual rules----------------------##
rule arcasHLA_extr_chr6:
     input:
         in_sortbamfile = "analysis/star/{sample}/{sample}.sorted.bam"
     output:
         chr6fastqfile1="analysis/neoantigen/{sample}/{sample}.sorted.extracted.1.fq.gz",
         chr6fastqfile2="analysis/neoantigen/{sample}/{sample}.sorted.extracted.2.fq.gz",
         outfile1 = "analysis/neoantigen/{sample}/{sample}.1.fq.gz",
         outfile2 = "analysis/neoantigen/{sample}/{sample}.2.fq.gz"
     threads: _neoantigen_individual_threads
     message: "Running ArcasHLA on {wildcards.sample}"
     log:
         "logs/neoantigen/{sample}.arcasHLA.log"
     params:
         arcasHLA_path=config["arcasHLA_path"],
         outpath = "analysis/neoantigen/{sample}",
         #outfile1 = "analysis/neoantigen/{sample}/{sample}.1.fq.gz",
         #outfile2 = "analysis/neoantigen/{sample}/{sample}.2.fq.gz"
     shell:
        """{params.arcasHLA_path}/arcasHLA extract {input.in_sortbamfile} --paired -t {threads} -v -o {params.outpath}"""
        """ && cp {output.chr6fastqfile1} {output.outfile1}"""
        """ && cp {output.chr6fastqfile2} {output.outfile2}"""

rule arcasHLA_genotype:
    input:
        fastq1 = "analysis/neoantigen/{sample}/{sample}.1.fq.gz",
        fastq2 = "analysis/neoantigen/{sample}/{sample}.2.fq.gz"
    output:
        "analysis/neoantigen/{sample}/{sample}.alignment.p",
        "analysis/neoantigen/{sample}/{sample}.genes.json",
        "analysis/neoantigen/{sample}/{sample}.genotype.json",
        "analysis/neoantigen/{sample}/{sample}.genotype.log"
    params:
        arcasHLA_path = config["arcasHLA_path"],
        outpath = "analysis/neoantigen/{sample}",
    shell:
        """{params.arcasHLA_path}/arcasHLA genotype {input.fastq1} {input.fastq2} -g A,B,C,DQA1,DQB1,DRB1 -t 16 -v -o {params.outpath} """


      
rule arcasHLA_relocate:
    input:
        "analysis/neoantigen/{sample}/{sample}.genotype.json",
    output:
        "analysis/neoantigen/merge/{sample}.genotype.json",
    params:
        outpath = "analysis/neoantigen/merge",
    shell:
        """cp {input} {params.outpath}"""


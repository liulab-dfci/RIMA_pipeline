#!/usr/bin/env python

#-------------------------------
# @author: Lin Yang; @ref: VIPER
# @email: jwang0611@gmail.com

_arcasHLA_threads=16 

def arcasHLA_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for sample in config["samples"]:
        ls.append("analysis/arcasHLA/%s/%s.alignment.p" % (sample,sample))
        ls.append("analysis/arcasHLA/%s/%s.genes.json" % (sample,sample))
        ls.append("analysis/arcasHLA/%s/%s.genotype.json" % (sample,sample))
        ls.append("analysis/arcasHLA/%s/%s.genotype.log" % (sample,sample))
        ls.append("analysis/arcasHLA/%s/%s.sorted.extracted.1.fq.gz" % (sample,sample))
        ls.append("analysis/arcasHLA/%s/%s.sorted.extracted.2.fq.gz" % (sample,sample))
        ls.append("analysis/arcasHLA/merge/%s.genotype.json" % (sample))
    ls.append("analysis/arcasHLA/merge/genotypes.tsv")
    ls.append("analysis/arcasHLA/merge/genotypes.p-group.tsv")
    ls.append("images/arcasHLA/hla_typing_frequency_plot.png")
    ls.append("files/multiqc/arcasHLA/hla_heatmap.txt")
    return ls

rule arcasHLA_all:
     input:
         arcasHLA_targets

rule arcasHLA_extr_chr6:
     input:
         in_sortbamfile = "analysis/star/{sample}/{sample}.sorted.bam"
     output:
         chr6fastqfile1="analysis/arcasHLA/{sample}/{sample}.sorted.extracted.1.fq.gz",
         chr6fastqfile2="analysis/arcasHLA/{sample}/{sample}.sorted.extracted.2.fq.gz"
     threads: _arcasHLA_threads
     message: "Running ArcasHLA on {wildcards.sample}"
     log:
         "logs/arcasHLA/{sample}.arcasHLA.log"
     params:
         arcasHLA_path=config["arcasHLA_path"],
         outpath = "analysis/arcasHLA/{sample}",
     shell:
        """{params.arcasHLA_path}/arcasHLA extract {input.in_sortbamfile} --paired -t {threads} -v -o {params.outpath}"""

rule arcasHLA_genotype:
    input:
        fastq1 = "analysis/arcasHLA/{sample}/{sample}.sorted.extracted.1.fq.gz",
        fastq2 = "analysis/arcasHLA/{sample}/{sample}.sorted.extracted.2.fq.gz"
    output:
        "analysis/arcasHLA/{sample}/{sample}.alignment.p",
        "analysis/arcasHLA/{sample}/{sample}.genes.json",
        "analysis/arcasHLA/{sample}/{sample}.genotype.json",
        "analysis/arcasHLA/{sample}/{sample}.genotype.log"
    params:
        arcasHLA_path=config["arcasHLA_path"],
        outpath = "analysis/arcasHLA/{sample}",
    shell:
        """{params.arcasHLA_path}/arcasHLA genotype {input.fastq1} {input.fastq2} -g A,B,C,DQA1,DQB1,DRB1 -t 16 -v -o {params.outpath}"""

rule arcasHLA_relocate:
    input:
        "analysis/arcasHLA/{sample}/{sample}.genotype.json",
    output:
        "analysis/arcasHLA/merge/{sample}.genotype.json",
    params:
        outpath = "analysis/arcasHLA/merge",
    shell:
        """cp {input} {params.outpath}"""

rule arcasHLA_merge:
    input:
        re = expand("analysis/arcasHLA/merge/{sample}.genotype.json", sample = config['samples'])
    output:
        "analysis/arcasHLA/merge/genotypes.tsv"
    params:
        arcasHLA_path=config["arcasHLA_path"],
        outpath = "analysis/arcasHLA/merge",
    shell:
        """{params.arcasHLA_path}/arcasHLA merge -i {params.outpath} -o {params.outpath}"""

rule arcasHLA_convert:
    input:
        res = "analysis/arcasHLA/merge/genotypes.tsv",
    output:
        pgroup = "analysis/arcasHLA/merge/genotypes.p-group.tsv",
    params:
        arcasHLA_path=config["arcasHLA_path"],
    shell:
        """{params.arcasHLA_path}/arcasHLA convert -r p-group {input.res} -o {output.pgroup}"""

rule arcasHLA_plot:
    input:
        res = "analysis/arcasHLA/merge/genotypes.p-group.tsv",
        expr = "analysis/combat/tpm_convertID.batch"
    output:
        png = "images/arcasHLA/hla_typing_frequency_plot.png",
        arcasHLA_table = "files/multiqc/arcasHLA/hla_heatmap.txt"
    params:
        meta = config["metasheet"],
        group = lambda wildcards: ','.join(str(i) for i in config["hla_annot_group"]),
        multiqc = "files/multiqc/arcasHLA/",
        outpath = "images/arcasHLA/",
        path="set +eu;source activate %s" % config['stat_root'],
    conda: "../../envs/stat_perl_r.yml"
    shell:
        "{params.path}; Rscript src/hla_plot.R --hla {input.res} --meta {params.meta} --expression {input.expr} --group {params.group} --outdir {params.outpath} --multiqc {params.multiqc}"

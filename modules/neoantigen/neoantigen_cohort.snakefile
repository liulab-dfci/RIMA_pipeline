#!/usr/bin/env python

design = config["design"]
batch = config["batch"]

#####################Neoantigen cohort module#####################################

def neoantigen_cohort_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    ls.append("analysis/neoantigen/merge/genotypes.tsv")
    ls.append("analysis/neoantigen/merge/genotypes.p-group.tsv")
    ls.append("analysis/neoantigen/%s_%s_hla_frequency_plot.png" % (design,batch))
    return ls

rule neoantigen_all:
     input:
         neoantigen_cohort_targets

##--------------------------arcasHLA cohort rules-----------------##
rule arcasHLA_merge:
    input:
        re = expand("analysis/neoantigen/merge/{sample}.genotype.json", sample = config['samples'])
    output:
        "analysis/neoantigen/merge/genotypes.tsv"
    params:
        arcasHLA_path=config["arcasHLA_path"],
        outpath = "analysis/neoantigen/merge",
    shell:
        """{params.arcasHLA_path}/arcasHLA merge -i {params.outpath} -o {params.outpath}"""

rule arcasHLA_convert:
    input:
        res = "analysis/neoantigen/merge/genotypes.tsv",
    output:
        pgroup = "analysis/neoantigen/merge/genotypes.p-group.tsv",
    params:
        arcasHLA_path=config["arcasHLA_path"],
    shell:
        """{params.arcasHLA_path}/arcasHLA convert -r p-group {input.res} -o {output.pgroup}"""
        
rule arcasHLA_plot:
    input:
        res = "analysis/neoantigen/merge/genotypes.p-group.tsv",
        expr = "analysis/batchremoval/tpm.genesymbol.batchremoved.csv"
    output:
        png = "analysis/neoantigen/{design}_{batch}_hla_frequency_plot.png"
    params:
        meta = config["metasheet"],
        design = config["design"],
        outpath = "analysis/neoantigen/{design}_{batch}_",
        path="set +eu;source activate %s" % config['stat_root'],
    conda: "../../envs/stat_perl_r.yml"
    shell:
        "{params.path}; Rscript src/neoantigen/hla_plot.R --hla {input.res} --meta {params.meta} --expression {input.expr} --design {params.design} --outdir {params.outpath}"

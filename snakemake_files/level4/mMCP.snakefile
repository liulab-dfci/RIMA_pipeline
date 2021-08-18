#!/usr/bin/env python

#-------------------------------
# @author: Jin Wang; @ref: VIPER
# @email: jwang0611@gmail.com


def mMCP_targets(wildcards):
    ls = []
    ls.append("files/mMCP/mMCP.txt")
    return ls

rule mMCP_all:
    input:
        mMCP_targets

rule mMCP_infiltration:
    input:
        "analysis/combat/tpm_convertID.txt"
    output:
        "files/mMCP/mMCP.txt"
    log:
        "logs/mMCP/mMCP.log"
    params:
        out_dir = "files/mMCP",
        path="set +eu;source activate %s" % config['mMCP_env']
    message: 
        "Running mMCP counter on the expression data"
    shell:
        "{params.path}; Rscript src/mMCP.R -e {input}  -o {params.out_dir}"




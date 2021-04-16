#!/usr/bin/env python

#-------------------------------Conda environment generation----------------------#

def environment_targets(wildcards):
    ls = []
    ls.append("files/environment/environment_log.txt")
    return ls


rule environment_all:
  input:
    environment_targets


rule environment_run:
 input:
    data = "static/environments/"
 output:
    out_file = "files/environment/environment_log.txt"
 log:
    "logs/environment/environments.txt"
 shell:
    "bash src/environment/create_environment.sh > {output.out_file}"

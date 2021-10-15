#!/usr/bin/env python
#----------------main snakefile to run RIMA-------------------------#
import os
import sys
import subprocess
import pandas as pd
import yaml

from string import Template
#from snakemake_files.scripts.utils import getTargetInfo

def getRuns(config):
    ret = {}
    #LEN: Weird, but using pandas to handle the comments in the file
    #KEY: need skipinitialspace to make it fault tolerant to spaces!
    metadata = pd.read_csv(config['metasheet'], index_col=0, sep=',', comment='#', skipinitialspace=True)
    f = metadata.to_csv().split() #make it resemble an actual file with lines
    #SKIP the hdr
    for l in f[1:]:
        tmp = l.strip().split(",")
        #print(tmp)
        ret[tmp[0]] = tmp[1:]
        #print(ret)
        config['runs'] = ret
    return config


def addCondaPaths_Config(config):
    """ADDS the python2 paths to config"""
    conda_root = subprocess.check_output('conda info --root',shell=True).decode('utf-8').strip()
    conda_path = os.path.join(conda_root, 'pkgs')
    current_path = os.getcwd()
    config['conda_root'] = conda_root
    config['stat_root'] = "%s/envs/stat_perl_r" % conda_root
    config['centrifuge_root']= "/%s/envs/centrifuge_env" % conda_root
    config['rseqc_root']= "/%s/envs/rseqc_env" % conda_root


def addExecPaths(config):
    conda_root = subprocess.check_output('conda info --root',shell=True).decode('utf-8').strip()
    conda_path = os.path.join(conda_root, 'pkgs')
    current_path = os.getcwd()
    #NEED the following when invoking python2 (to set proper PYTHONPATH)
    if not "trust4_path" in config or not config["trust4_path"]:
        config["trust4_path"] = os.path.join(current_path,'TRUST4')
    return config

def loadRef(config):
    f = open(config['ref'])
    ref_info = yaml.safe_load(f)
    f.close()
    #print(ref_info[config['assembly']])
    for (k,v) in ref_info[config['assembly']].items():
    #NO CLOBBERING what is user-defined!
        if k not in config:
            config[k] = v

def load_config(config_file):
    #load the main config file including parameters which are not change a lot
    with open(config_file, 'r') as stream:
        try:
            return yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

def load_execution(execution_file):
    #load the main config file including parameters which are not change a lot
    with open(execution_file, 'r') as stream:
        try:
            return yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

#---------  CONFIG set up  ---------------
config = load_config('config.yaml')
execution = load_execution('execution.yaml')
addCondaPaths_Config(config)
loadRef(config)
addExecPaths(config)




#------------------------------------------------------------------------------
# TARGETS
#------------------------------------------------------------------------------
def all_targets(wildcards):
    ls = []
    if execution["preprocess_individual"]:
        ls.extend(preprocess_individual_targets(wildcards))
    if execution["preprocess_cohort"]:
        ls.extend(preprocess_cohort_targets(wildcards))
    if execution["differential_expression_cohort"]:
        ls.extend(diffexpr_targets(wildcards))
    if execution["immune_infiltration_cohort"]:
        ls.extend(immune_infiltration_targets(wildcards))
    if execution["immune_repertoire_individual"]:
        ls.extend(immune_repertoire_individual_targets(wildcards))
    if execution["microbiome_individual"]:
        ls.extend(microbiome_individual_targets(wildcards))
    if execution["microbiome_cohort"]:
        ls.extend(microbiome_cohort_targets(wildcards))    
    return ls


rule target:
    input:
        all_targets
    message: "Compiling all outputs"


if execution["preprocess_individual"]:
    include: "./modules/preprocess/preprocess_individual.snakefile"
if execution["preprocess_cohort"]:
    include: "./modules/preprocess/preprocess_cohort.snakefile"
if execution["differential_expression_cohort"]:
    include: "./modules/differential_expression/differential_expression_cohort.snakefile"
if execution["immune_infiltration_cohort"]:
    include: "./modules/immune_infiltration/immune_infiltration_cohort.snakefile"
if execution["immune_repertoire_individual"]:
    include: "./modules/immune_repertoire/immune_repertoire_individual.snakefile"
if execution["microbiome_individual"]:
    include: "./modules/microbiome/microbiome_individual.snakefile"
if execution["microbiome_cohort"]:
    include: "./modules/microbiome/microbiome_cohort.snakefile"

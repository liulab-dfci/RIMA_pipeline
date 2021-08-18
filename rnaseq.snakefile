#!/usr/bin/env python

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

# def addCondaPaths_Config(config):
#     """ADDS the python2 paths to config"""
#     conda_root = subprocess.check_output('conda info --root',shell=True).decode('utf-8').strip()
#     config['conda_root'] = conda_root
#     config['wes_root'] = "%s/envs/wes" % conda_root

def addCondaPaths_Config(config):
    """ADDS the python2 paths to config"""
    conda_root = subprocess.check_output('conda info --root',shell=True).decode('utf-8').strip()
    conda_path = os.path.join(conda_root, 'pkgs')
    current_path = os.getcwd()
    config['conda_root'] = conda_root
    config['rnaseq_root'] = "%s/envs/rna" % conda_root
    config['stat_root'] = "%s/envs/stat_perl_r" % conda_root
    config['centrifuge_root']= "/%s/envs/centrifuge_env" % conda_root
    config['rseqc_root']= "%s/envs/rseqc_env" % conda_root
    config['gatk4_root']= "/%s/envs/gatk4_env" % conda_root
    config['msisensor_root']= "/%s/envs/msisensor_env" % conda_root
    config['vep_root']= "/%s/envs/vep_env" % conda_root
    config['varscan_root']= "/%s/envs/varscan_env" % conda_root
    config['immnue_decov']= "%s/envs/imm_env" % conda_root
    config['trust4_env']= "%s/envs/trust4_env" % conda_root
    config['mMCP_env']= "%s/envs/mMCP_env" % conda_root
    config['pvacseq'] = "%s/envs/pvacseq" % conda_root

    
def addExecPaths(config):
    conda_root = subprocess.check_output('conda info --root',shell=True).decode('utf-8').strip()
    conda_path = os.path.join(conda_root, 'pkgs')
    current_path = os.getcwd()
    TRUST4_path = "/liulab/RIMA_kraken"
    arcasHLA_path = "/liulab/linyang"
    #NEED the following when invoking python2 (to set proper PYTHONPATH)
    if not "pvacseq_path" in config or not config["pvacseq_path"]:
        config["pvacseq_path"] = os.path.join(conda_root, 'envs', 'pvacseq', 'bin')
    if not "trust4_path" in config or not config["trust4_path"]:
        config["trust4_path"] = os.path.join(TRUST4_path,'TRUST4')
    if not "arcasHLA_path" in config or not config["arcasHLA_path"]:
        config["arcasHLA_path"] = os.path.join(arcasHLA_path,'arcasHLA')
    if not "msisensor2_path" in config or not config["msisensor2_path"]:
        config["msisensor2_path"] = os.path.join(TRUST4_path,'msisensor2')
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
#config = getRuns(config)
loadRef(config)
addExecPaths(config)




#------------------------------------------------------------------------------
# TARGETS
#------------------------------------------------------------------------------
def all_targets(wildcards):
    ls = []
    #------------------------------------------------------------------------
    ##Level1
    #------------------------------------------------------------------------
    if execution["star"]:
        ls.extend(align_targets(wildcards))
        ls.extend(star_report_tagets(wildcards))
    if execution["salmon"]:
        ls.extend(salmon_targets(wildcards))
    if execution["rseqc"]:
        ls.extend(rseqc_targets(wildcards))
    #-----------------------------------------------------------------------
    ##Level2
    #----------------------------------------------------------------------
    if execution["batch_removal"]:
        ls.extend(combat_targets(wildcards))
    if execution["deseq2"]:
        ls.extend(deseq2_targets(wildcards))
    if execution["gsea"]:
        ls.extend(gsea_targets(wildcards))
    if execution["ssgsea"]:
        ls.extend(ssgsea_targets(wildcards))
    #--------------------------------------------------------------------    
    ##Level3
    #--------------------------------------------------------------------
 #   if execution["fusion"]:
 #       ls.extend(fusion_targets(wildcards))
    if execution["microbiota"]:
        ls.extend(microbe_targets(wildcards))
    if execution["trust4"]:
        ls.extend(trust4_targets(wildcards))
    if execution["arcasHLA"]:
        ls.extend(arcasHLA_targets(wildcards))
    #---------------------------------------------------------------------    
    ##Level4
    #---------------------------------------------------------------------
#    if execution["tide"]:
#        ls.extend(tide_targets(wildcards))
    if execution["immunedeconv"]:
        ls.extend(immunedeconv_targets(wildcards))
    if execution["mMCP"]:
        ls.extend(mMCP_targets(wildcards))
    #-------------------------------------------------------------------
    ##Level5
    #-------------------------------------------------------------------
    if execution["msisensor2"]:
        ls.extend(msisensor_targets(wildcards))
    if execution["report"]:
        ls.extend(report_targets(wildcards))
    if execution["pvacseq"]:
        ls.extend(neoantigen_targets(wildcards))

    return ls

rule target:
    input: 
        all_targets
    message: "Compiling all output"

 ##Level1   
if execution["star"]:
    include: "./snakemake_files/level1/star_final.snakefile"
if execution["salmon"]:
    include: "./snakemake_files/level1/salmon_final.snakefile"
if execution["rseqc"]:
    include: "./snakemake_files/level1/rseqc_final.snakefile"
##Level2
if execution["batch_removal"]:
    include: "./snakemake_files/level2/batch_removal.snakefile"
if execution["deseq2"]:
    include: "./snakemake_files/level2/deseq2.modified.snakefile"
if execution["gsea"]:
    include: "./snakemake_files/level2/gsea.modified.snakefile"
if execution["ssgsea"]:
    include: "./snakemake_files/level2/ssgsea_final.snakefile"
##Level3
#if execution["fusion"]:
#    include: "./snakemake_files/level3_individual.snakefile"
if execution["microbiota"]:
    include: "./snakemake_files/level3/microbiota_final.snakefile"
if execution["trust4"]:
    include: "./snakemake_files/level3/trust4_final.snakefile"
if execution["arcasHLA"]:
    include: "./snakemake_files/level3/arcasHLA.snakefile"
##Level4
#if execution["tide"]:
#    include: "./snakemake_files/level4/tide.snakefile"
if execution["immunedeconv"]:
    include: "./snakemake_files/level4/immunedeconv.snakefile"
if execution["mMCP"]:
    include: "./snakemake_files/level4/mMCP.snakefile"
##Level5
if execution["msisensor2"]:
    include: "./snakemake_files/level5/msisensor.snakefile"
if execution["report"]:
    include: "./snakemake_files/report.snakefile"
if execution["pvacseq"]:
    include: "./snakemake_files/level5/pvacseq_final.snakefile"

    














#!/usr/bin/env python

#-------------------------------

metadata = pd.read_csv(config["metasheet"], index_col=0, sep=',')
options = [config["Treatment"],config["Control"]]
design = config["design"]
treatment = config["Treatment"]
control = config["Control"] 
pretreat = config["pre_treated"]
cancer = config["cancer_type"]

if cancer in ['NSCLC','Melanoma']:
	cancer_type = cancer
else:
	cancer_type = 'Other'
	
def getsampleIDs(meta):
	return meta[meta[design].isin(options)].index


def immune_response_cohort_targets(wildcards):
    ls = []
    ls.append("analysis/TIDE/%s_tide_input.txt" % design )
    ls.append("analysis/TIDE/%s_tide_output.txt" % design )
    ls.append("analysis/TIDE/%s_TIDE-TCGA_mqc.png" % design)
    ls.append("analysis/msisensor/%s_msi_score.txt" % design)
    ls.append("analysis/TIDE/%s_comparison.png" % design)

    return ls

rule immune_response_cohort_all:
    input:
        immune_response_cohort_targets

#------------------------TIDEpy rules-----------------------------#
rule immune_response_input:
    input:
        "analysis/batchremoval/tpm.genesymbol.batchremoved.csv"
    output:
        "analysis/TIDE/{design}_tide_input.txt"
    message:
        "Prepare TIDE input"
    benchmark:
        "benchmarks/immune_response/{design}_tide_input.benchmark"
    log:
        "logs/immune_response/{design}_tide_input.log"
    params:
        outdir = "analysis/TIDE/",
        cancer_type = cancer_type,
        design = design,
        pretreat = pretreat,
        path="set +eu;source activate %s" % config['stat_root'],
    shell:
    	"{params.path}; Rscript src/immune_response/tide_run.R --input {input} --design {params.design} \
    	--cancer {params.cancer_type} --treated {params.pretreat} --outdir {params.outdir}"
    	

rule immune_response_output:
    input:
        "analysis/TIDE/{design}_tide_input.txt"
    output:
        "analysis/TIDE/{design}_tide_output.txt"
    message:
        "Running TIDEpy"
    benchmark:
        "benchmarks/immune_response/{design}_tide_score.benchmark"
    log:
        "logs/immune_response/{design}_tide_score.log"
    params:
        outdir = "analysis/TIDE/",
        cancer_type = cancer_type,
        design = design,
        path="set +eu;source activate %s" % config['stat_root'],
    run:
        if pretreat == 'True':
            shell("""tidepy -c {params.cancer_type} --pretreat -o {output} {input} """)
        else:
            shell("""tidepy -c {params.cancer_type} -o {output} {input} """)
    	
    	
rule immune_response_plot:
    input:
        score = "analysis/TIDE/{design}_tide_output.txt",
        expr ="analysis/TIDE/{design}_tide_input.txt"
    output:
        "analysis/TIDE/{design}_TIDE-TCGA_mqc.png"
    message:
        "plot on tide score"
    benchmark:
        "benchmarks/immune_response/{design}_tide_plot.benchmark"
    log:
        "logs/immune_response/{design}_tide_plot.log"
    params:
        cancer = config["cancer_type"], 
        outpath = "analysis/TIDE/",
        design = design,
        path="set +eu;source activate %s" % config['stat_root']
    conda:
        "../envs/stat_perl_r.yml"
    shell:
        "{params.path}; Rscript src/immune_response/tide_plot.R --input {input.score} \
        --design {params.design} --expression {input.expr} --cc {params.cancer} --outdir {params.outpath}"
        

rule merge_msisensor:
    input:
        files = expand("analysis/msisensor/{sample}/{sample}_msisensor", sample = getsampleIDs(metadata) )
    output:
        "analysis/msisensor/{design}_msi_score.txt"
    message:
    	"Merging msisensor scores "
    benchmark:
    	"benchmarks/immune_response/{design}_msi_plot.benchmark"
    log:
    	"analysis/variant/{design}_variant_missensor.log"
    conda: 
    	"../envs/stat_perl_r.yml"
    params:
    	outpath = "analysis/msisensor/",
    	filelist = lambda wildcards, input: ','.join(str(i) for i in list({input.files})[0]),
    	design = design,
    	meta = config["metasheet"],
    	path="set +eu;source activate %s" % config['stat_root']
    shell:
    	"{params.path}; Rscript src/immune_response/merge_msi.R --input {params.filelist} \
    	--meta {params.meta} --outdir {params.outpath} --condition {params.design}"


rule immune_comparison_plot:
    input:
        msi = "analysis/msisensor/{design}_msi_score.txt",
        tide = "analysis/TIDE/{design}_tide_output.txt"
    output:
        "analysis/TIDE/{design}_comparison.png"
    log:
        "logs/TIDE/{design}_compare.log"
    message:
        "Running immune response ploting"
    benchmark:
        "benchmarks/immune_response/{design}_compare.benchmark"
    conda: "../envs/stat_perl_r.yml"
    params:
        outpath = "analysis/TIDE/",
        condition = design,
        treatment = treatment,
        control = control,
        meta = config["metasheet"],
        path="set +eu;source activate %s" % config['stat_root'],
    shell:
        "{params.path};Rscript src/immune_response/response_plot.R --msiscore {input.msi} --tidescore {input.tide} \
        --meta {params.meta} --outdir {params.outpath} --condition {params.condition} --treatment {params.treatment} --control {params.control}"

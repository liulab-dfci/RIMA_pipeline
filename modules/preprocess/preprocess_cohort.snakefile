#!/usr/bin/env python

#-------------------------------Preprocessing cohort module------------------------#
_preprocess_threads = 8

#-------------------Preprocess cohort targets----------------------#
import pandas as pd

design = config["design"]
covariates = config["batch"]


def merge_sep_inputs(inputs):
    inputs_format = ' -f '.join(str(i) for i in list(inputs)[0])
    return inputs_format

def merge_sal_inputs(inputs):
    inputs_format = ','.join(str(i) for i in list(inputs)[0])
    return inputs_format


def preprocess_cohort_targets(wildcards):
    ls = []
    ls.append("analysis/star/STAR_Align_Report.csv" )
    ls.append("analysis/rseqc/gene_body_cvg/geneBodyCoverage.r")
    ls.append("analysis/rseqc/gene_body_cvg/geneBodyCoverage.curves.pdf")
    ls.append("analysis/rseqc/tin_score/tin_score_summary.txt")
    ls.append("analysis/rseqc/read_distrib/read_distrib.matrix.tab")
    
    ls.append("analysis/salmon/tpm.genesymbol.csv")
    ls.append("analysis/batchremoval/%s_%s_tpm.genesymbol.csv" % (design,covariates))
    ls.append("analysis/batchremoval/%s_%s_tpm.genesymbol.batchremoved.csv" % (design,covariates))
    
    ls.append("analysis/batchremoval/%s_%s_pca_plot_before.pdf" % (design,covariates))
    ls.append("analysis/batchremoval/%s_%s_pca_plot_after.pdf" % (design,covariates))
    return ls 

rule preprocess_cohort_all:
    input:
      preprocess_cohort_targets

#--------------------STAR_cohort_csv---------------------------#
rule STAR_matrix:
    input:
      star_log_files=expand( "analysis/star/{sample}/{sample}.Log.final.out", sample=config["samples"] )
    output:
      csv="analysis/star/STAR_Align_Report.csv"
    message:
      "Generating STAR report"
    benchmark:
      "benchmarks/star/generate_STAR_report.benchmark"
    conda: "../envs/stat_perl_r.yml"
    params:
      log_files = lambda wildcards, input: merge_sep_inputs({input.star_log_files}),
      path="set +eu;source activate %s" % config['stat_root'],
    shell:
      """{params.path}; perl src/preprocess/STAR_reports.pl -f {params.log_files} 1>{output.csv} """

#--------------------------RSeQC cohort----------------#
rule tin_summary:
    input:
      expand("analysis/rseqc/tin_score/{sample}/{sample}.summary.txt", sample=config['samples'])
    output:
      score = "analysis/rseqc/tin_score/tin_score_summary.txt",
    message:
      "plotting TIN score summary"
    log:
      "logs/rseqc/tin_score/tin_score_summary.log"
    benchmark:
      "benchmarks/rseqc/tin_score/tin_score_summary.benchmark"
    params:
      path="set +eu;source activate %s" % config['stat_root'],
    conda: "../envs/stat_perl_r.yml"
    shell:
      """cat {input} | sed '1 !{{/Bam_file/d;}}' >{output.score}"""

rule read_distrib_qc_matrix:
    input:
      read_distrib_files=expand( "analysis/rseqc/read_distrib/{sample}/{sample}.txt", sample = config['samples'])
    output:
      matrix="analysis/rseqc/read_distrib/read_distrib.matrix.tab",
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
      "perl src/preprocess/read_distrib_matrix.pl -f {params.file_list_with_flag} 1>{output.matrix} "

rule plot_gene_body_cvg:
    input:
      samples_list=expand("analysis/rseqc/gene_body_cvg/{sample}/{sample}.geneBodyCoverage.r", sample=config["samples"] )
    output:
      rscript="analysis/rseqc/gene_body_cvg/geneBodyCoverage.r",
      png_curves="analysis/rseqc/gene_body_cvg/geneBodyCoverage.curves.pdf"
    message: "Plotting gene body coverage"
    benchmark:
      "benchmarks/rseqc/gene_body_cvg/plot_gene_body_cvg.benchmark"
    params:
      path="set +eu;source activate %s" % config['stat_root'],
    conda: "../envs/stat_perl_r.yml"
    shell:
      "{params.path};perl src/preprocess/plot_gene_body_cvg.pl --rfile {output.rscript} --curves_png {output.png_curves}"
      
      " {input.samples_list} && {params.path}; Rscript {output.rscript}"
      


#-------------Salmon cohort csv----------------------#
rule salmon_matrix:
    input:
      salmon_tpm_files = expand("analysis/salmon/{sample}/{sample}.quant.sf", sample = config['samples']),
      metasheet = config["metasheet"]
    output:
      "analysis/salmon/tpm.genesymbol.csv"
    benchmark:
      "benchmarks/salmon/salmon_gene_matrix.benchmark"
    params:
      args = lambda wildcards, input: merge_sal_inputs({input.salmon_tpm_files}),
      tx2gene = 'ref_files/tximport/tx2gene.csv',
      outpath = 'analysis/salmon/',
      path = "set +eu;source activate %s" % config['stat_root']
      
    message: "Merge Salmon gene quantification together for all samples "
    shell:
      "{params.path}; Rscript src/preprocess/merge_tpm.R --input {params.args} --meta {input.metasheet} \
      --type salmon --tx2gene {params.tx2gene} --outpath {params.outpath}"
      
      

#--------------------Batch removal---------------------#
rule batch_removal:
    input:
        "analysis/salmon/tpm.genesymbol.csv"
    output:
        after = "analysis/batchremoval/{design}_{covariates}_tpm.genesymbol.batchremoved.csv",
        before = "analysis/batchremoval/{design}_{covariates}_tpm.genesymbol.csv",
    message:
        "Running batch removal using limma method"
    benchmark:
        "benchmarks/batchremoval/{design}_{covariates}_tpm_limma.benchmark" 
    params:
        covariates = config["batch"],
        design = config["design"],
        path="set +eu;source activate %s" % config['stat_root'],
        meta = config["metasheet"],
        rename = "analysis/batchremoval/tpm.genesymbol.batchremoved.csv"
    log:
        "logs/batchremoval/{design}_{covariates}_batch_removal.log"
    conda: "../envs/stat_perl_r.yml"
    shell:
        "{params.path}; Rscript src/preprocess/batch_removal.R -e {input} -c {params.covariates} \
        -d {params.design} -m {params.meta} -b {output.before} -a {output.after} \
        && cp {output.after} {params.rename}"


        
rule pca_sample_clustering:
    input:
        after_batch = "analysis/batchremoval/{design}_{covariates}_tpm.genesymbol.batchremoved.csv",
        before_batch = "analysis/batchremoval/{design}_{covariates}_tpm.genesymbol.csv"
    output:
        before_pca = "analysis/batchremoval/{design}_{covariates}_pca_plot_before.pdf",
        after_pca = "analysis/batchremoval/{design}_{covariates}_pca_plot_after.pdf",
    message:
        "Running PCA for sample clustering"
    benchmark:
        "benchmarks/batchremoval/{design}_{covariates}_pca.benchmark"
    params:
        meta_info = config["metasheet"],
        path="set +eu;source activate %s" % config['stat_root'],
        covariates = config["batch"],
        design = config["design"]
    conda: "../envs/stat_perl_r.yml"
    shell:
        "{params.path}; Rscript src/preprocess/pca.R -b {input.before_batch} -a {input.after_batch} -m {params.meta_info}  -c {params.covariates} \
                -g {params.design} -i {output.before_pca} -j {output.after_pca}"



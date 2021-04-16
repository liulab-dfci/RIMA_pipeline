#!/usr/bin/env python

#-------------------------------Preprocessing cohort module------------------------#
_preprocess_threads = 8

#-------------------Preprocess cohort targets----------------------#
import pandas as pd
metadata = pd.read_csv(config["metasheet"], index_col=0, sep=',')


##extract gene reads count column from star output
if config["library_type"] == "fr-firststrand":
    count_col = 3
elif config["library_type"] == "fr-secondstrand":
    count_col = 2
else:
    count_col = 1



def merge_sep_inputs(inputs):
    inputs_format = ' -f '.join(str(i) for i in list(inputs)[0])
    return inputs_format

def preprocess_cohort_targets(wildcards):
    ls = []
    ls.append("files/star/STAR_Align_Report.csv" )
    ls.append("analysis/star/STAR_Gene_Counts.csv")
    ls.append("analysis/rseqc/gene_body_cvg/geneBodyCoverage.r")
    ls.append("files/rseqc/gene_body_cvg/geneBodyCoverage.curves.png")
    ls.append("files/rseqc/tin_score/tin_score_summary.txt")
    ls.append("analysis/rseqc/read_distrib/read_distrib.matrix.tab")
    ls.append("analysis/salmon/salmon_tpm_matrix.csv")
    ls.append("analysis/batchremoval/tpm_matrix.batch")
    ls.append("analysis/batchremoval/tpm_convertID.txt")
    ls.append("analysis/batchremoval/tpm_convertID.batch")
    ls.append("analysis/batchremoval/gencode_tx2_ensemble_gene_symbol.csv")
    ls.append("files/batchremoval/pca_plot_before.png")
    ls.append("files/batchremoval/pca_plot_after.png")
    ls.append("analysis/batchremoval/tpm_convertID_batch_Entrez.txt")
    ls.append("analysis/batchremoval/tide/tpm_convertID_batch_Entrez.txt")
    ls.append("analysis/batchremoval/tide/tide_meta.txt")
    return ls

rule preprocess_cohort_all:
    input:
      preprocess_cohort_targets

#--------------------STAR_cohort_csv---------------------------#
rule STAR_matrix:
    input:
      star_log_files=expand( "analysis/star/{sample}/{sample}.Log.final.out", sample=config["samples"] ),
      star_gene_count_files=expand( "analysis/star/{sample}/{sample}.counts.tab", sample=config["samples"] )
    output:
      csv="files/star/STAR_Align_Report.csv",
      gene_counts="analysis/star/STAR_Gene_Counts.csv"
    message:
      "Generating STAR report"
    benchmark:
      "benchmarks/star/generate_STAR_report.benchmark"
    conda: "../envs/stat_perl_r.yml"
    params:
      extract_col = count_col,
      log_files = lambda wildcards, input: merge_sep_inputs({input.star_log_files}),
      count_files = lambda wildcards, input: merge_sep_inputs({input.star_gene_count_files}),
      path="set +eu;source activate %s" % config['stat_root'],
    shell:
      """{params.path}; perl src/preprocess/STAR_reports.pl -f {params.log_files} 1>{output.csv} && """
      """{params.path}; perl src/preprocess/raw_count_fpkm_tpm_matrix.pl --column {params.extract_col} -f {params.count_files} 1>{output.gene_counts}"""

#--------------------------RSeQC cohort----------------#
rule tin_summary:
    input:
      expand("analysis/rseqc/tin_score/{sample}/{sample}.summary.txt", sample=config['samples'])
    output:
      score = "files/rseqc/tin_score/tin_score_summary.txt",
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
      samples_list=expand("files/rseqc/gene_body_cvg/{sample}/{sample}.geneBodyCoverage.r", sample=config["samples"] )
    output:
      rscript="analysis/rseqc/gene_body_cvg/geneBodyCoverage.r",
      png_curves="files/rseqc/gene_body_cvg/geneBodyCoverage.curves.png"
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
      "analysis/salmon/salmon_tpm_matrix.csv"
    benchmark:
      "benchmarks/salmon/salmon_gene_matrix.benchmark"
    params:
      args = lambda wildcards, input: merge_sep_inputs({input.salmon_tpm_files})
    message: "Merge Salmon gene quantification together for all samples "
    shell:
      "perl src/preprocess/raw_count_fpkm_tpm_matrix.pl --metasheet {input.metasheet} -s --header -f {params.args} 1>{output}"

#--------------------Batch removal---------------------#

rule batch_removal:
    input:
        "analysis/salmon/salmon_tpm_matrix.csv"
    output:
        combat_expr = "analysis/batchremoval/tpm_matrix.batch"
    message:
        "Running batch removal using limma method"
    benchmark:
        "benchmarks/batchremoval/tpm_combat.benchmark"
    params:
        batch_file = config["metasheet"],
        covariates = lambda wildcards: ','.join(str(i) for i in config["batch_covariates"]),
        prefix = "analysis/batchremoval/tpm_matrix",
        path="set +eu;source activate %s" % config['stat_root'],
        #method = "combat"
    log:
        "logs/batchremoval/batch_removal.log"
    conda: "../envs/stat_perl_r.yml"
    shell:
        "{params.path}; Rscript src/preprocess/batch_removal.R -e {input} -b {params.batch_file} -c {params.covariates} -o {params.prefix}"

##changre 10 to 8 in 2nd cmd for v27        
rule id_convert:
    input:
        batch = "analysis/batchremoval/tpm_matrix.batch",
        raw = "analysis/salmon/salmon_tpm_matrix.csv"
    output:
        gene_convertID = "analysis/batchremoval/tpm_convertID.txt",
        after_batch = "analysis/batchremoval/tpm_convertID.batch",
        annotation = "analysis/batchremoval/gencode_tx2_ensemble_gene_symbol.csv"
    message:
        "Convert Ensemble ID or transcript ID to Gene Symbol from salmon matrix"
    log:
        "logs/batchremoval/convertID.log"
    benchmark:
        "benchmarks/batchremoval/convertID.benchmark"
    params:
        gene_prefix = "analysis/batchremoval/tpm_convertID.txt",
        tool_type = "salmon"
    shell:
        """zless {config[gtf_path]} | grep -v "#" | awk '$3=="transcript"' | cut -f9 | tr -s ";" " \
" | awk '{{print$4"\t"$2"\t"$8}}' | sort | uniq | sed 's/\"//g' > {output.annotation} """
        """ && python src/preprocess/ConvertID.py -m {input.batch} -t {params.tool_type} -r {output.annotation} -p {\
output.after_batch} """
        """ && python src/preprocess/ConvertID.py -m {input.raw} -t {params.tool_type} -r {output.annotation} -p {\
params.gene_prefix} """

        
rule pca_sample_clustering:
    input:
        before_batch = "analysis/batchremoval/tpm_convertID.txt",
        after_batch = "analysis/batchremoval/tpm_convertID.batch"
    output:
        "files/batchremoval/pca_plot_before.png",
        "files/batchremoval/pca_plot_after.png",
    message:
        "Running PCA for sample clustering"
    benchmark:
        "benchmarks/batchremoval/pca.benchmark"
    params:
        meta_info = config["metasheet"],
        out_path = "files/batchremoval/",
        path="set +eu;source activate %s" % config['stat_root'],
        covariates = lambda wildcards: ','.join(str(i) for i in config["batch_covariates"])
    conda: "../envs/stat_perl_r.yml"
    shell:
        "{params.path}; Rscript src/preprocess/pca.R -b {input.before_batch} -a {input.after_batch} -m {params.meta_info} -o {params.out_path} -c {params.covariates}"


#----------------separate into pre and post--------------------#

rule entrez_convert:
    input:
        "analysis/batchremoval/tpm_convertID.batch"
    output:
        entrez_file = "analysis/batchremoval/tpm_convertID_batch_Entrez.txt",
        tmp = "analysis/batchremoval/tide/tpm_convertID_batch_Entrez.txt"
    message:
        "Running ID entrez_convert"
    benchmark:
        "benchmarks/batchremoval/entrez_convert.benchmark"
    params:
        prefix = "analysis/batchremoval/tpm_convertID_batch_Entrez.txt",
        tool_type = "tide"
    shell:
        "python src/preprocess/ConvertID.py -m {input} -t {params.tool_type} -r static/msi_est/name.map -p {params.prefix}"
        " && tr ',' '\t' < {output.entrez_file} > {output.tmp}"

rule split_tpm_file:
   input:
        expr= "analysis/batchremoval/tide/tpm_convertID_batch_Entrez.txt"
   output:
       "analysis/batchremoval/tide/tide_meta.txt"
   params:
        meta_info = config["metasheet"],
        path = "set +eu;source activate %s" % config['stat_root'],
        outpath = "analysis/batchremoval/",
        #outpath2 = "analysis/batchremoval/post/"
   message:
        "Running preprocessing for pre and post separation"
   benchmark:
        "benchmarks/batchremoval/preprocessbatchfile.benchmark"
   conda:
        "../envs/stat_perl_r.yml"
   shell:
        "{params.path}; Rscript src/preprocess/new_batchparser.R -m {params.meta_info} -e {input.expr} --outdir {params.outpath}"



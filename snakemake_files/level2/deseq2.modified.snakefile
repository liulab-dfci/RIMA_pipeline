#File Edit Options Buffers Tools Python Help
#!/usr/bin/env python

#-------------------------------DESeq2 for differential gene expression-----------------------#
import itertools
import pandas as pd

data_type = config["assembly"]

def DiffInput(data_type):
        return expand("analysis/salmon/{sample}/{sample}.quant.sf",sample=config["samples"])

def DiffReference(data_type):
    if data_type == "mm10":
        reference = "static/deseq2/tx2gene_mouse.csv"
    if data_type == "hg38":
        reference = "static/deseq2/tx2gene.csv"
    return reference


def deseq2_targets(wildcards):
    ls=[]
    for design in config["designs"]:
        meta = config["metasheet"]
        design_file = "./" + meta
        design_meta = pd.read_csv(design_file, index_col=0, sep=',')
        design_meta = design_meta.rename(columns = {design: 'Condition'})
        comps = [i for i in list(set(design_meta["Condition"].tolist()))]
        combinations = list(itertools.combinations(comps,2))
        if len(comps) >1:
            if config["comparison"] == "between":
                compare_list = combinations
            if config["comparison"] == "loop":
                compare_list = comps
            for cp in compare_list:
                if config["comparison"] == "between":
                    cp = sorted(list(cp))
                    cp_list = cp[0]+"_VS_"+cp[1]
                else:
                    cp_list = cp+"_VS_others"
                ls.append("images/deseq2/%s/%s_%s_diff_volcano_plot.png" % (design,design,cp_list))
                ls.append("files/deseq2/%s/%s_%s_DESeq2_ConvertID.txt" % (design,design,cp_list))
                ls.append("analysis/deseq2/%s/%s_%s_DESeq2_raw.txt" % (design,design,cp_list))
                ls.append("files/multiqc/DESeq2/%s_%s_diff_volcano_plot.png" % (design,cp_list))
                ls.append("files/multiqc/DESeq2/%s_%s_DESeq2_sub.txt" % (design,cp_list))


    return ls


rule deseq2_differential_genes:
    input:
        files = DiffInput(data_type),
    output:
        deseq2_res = "files/deseq2/{design}/{design}_{compare}_DESeq2_ConvertID.txt",
        deseq2_raw = "analysis/deseq2/{design}/{design}_{compare}_DESeq2_raw.txt",
        multiqc = "files/multiqc/DESeq2/{design}_{compare}_DESeq2_sub.txt"
    log:
        "logs/deseq2/{design}/{design}_{compare}.deseq2.log"
    params:
        comps = config["comparison"],
        filelist = lambda wildcards, input: ','.join(str(i) for i in list({input.files})[0]),
        batch = config["batch"],
        out_path = "analysis/deseq2/{design}/{design}",
        file_path = "files/deseq2/{design}/",
        tx_annot = DiffReference(data_type),
        condition = config["designs"],
        meta = config["metasheet"],
        multiqc = "files/multiqc/DESeq2/{design}",
        source = "salmon",
        path="set +eu;source activate %s" % config['stat_root']
    message:
        "Running DESeq2 on the samples"
    benchmark:
        "benchmarks/deseq2/{design}/{design}_{compare}.deseq2.benchmark"
    conda: "../envs/stat_perl_r.yml"
    shell:
        "{params.path}; Rscript src/DESeq2.R --input {params.filelist} --batch {params.batch} --typ\
e {params.source} --comparison {params.comps} --meta {params.meta} --tx2gene {params.tx_annot} -\
-condition {params.condition} --outpath {params.out_path} --multiqc {params.multiqc}"
        " && mv {params.out_path}*DESeq2_ConvertID.txt {params.file_path}"


rule volcano_plot:
    input:
        "files/deseq2/{design}/{design}_{compare}_DESeq2_ConvertID.txt"
    output:
       plot =  "images/deseq2/{design}/{design}_{compare}_diff_volcano_plot.png",
       multiqc_plot = "files/multiqc/DESeq2/{design}_{compare}_diff_volcano_plot.png"
    params:
        path="set +eu;source activate %s" % config['stat_root'],
        multiqc_folder = "files/multiqc/DESeq2/"
    log:
        "logs/deseq2/{design}/{design}_{compare}.deseq2.volcano.log"
    message:
        "Running Volcano plot on the DESeq2 result"
    benchmark:
        "benchmarks/deseq2/{design}/{design}_{compare}.deseq2.volcano.benchmark"
    conda: "../envs/stat_perl_r.yml"
    shell:
        "{params.path}; Rscript src/volcano_plot.R --deseq2_mat {input} --outdir {output.plot}"
        " && cp {output.plot} {params.multiqc_folder}"
        " && cp {output.multiqc_plot} {params.multiqc_folder}DESeq2-Volcano-Plot_mqc.png"


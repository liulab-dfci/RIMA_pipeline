# How to Run RIMA on Kraken

## RIMA Workflow 
![](https://github.com/Lindky/RIMA_Kraken/blob/master/files/multiqc/Pipeline_Workflow_mqc.png)

### The full tutorial 
https://liulab-dfci.github.io/RIMA/


**Note:** The Kraken users can skip the RIMA installation step and downloading the reference. For the details, please follow this README.


## Available Tools Checklist
| **Methods** | **Description** | **Available models**|
| :---: | :---: | :---: |
| STAR | Spliced Transcripts Alignment to a Reference | Human, Mouse |
| Salmon | Gene Quantification | Human, Mouse |
| RSeQC | High Throughput Sequence Data Evaluation | Human, Mouse |
| batch_removal| Remove Batch Effects Using Limma | Human, Mouse |
|  | **---DIFFERENTIAL EXPRESSION---** |  |
| DESeq2 | Gene Differential Expression Analysis | Human, Mouse |
| GSEA | Gene Set Enrichment Analysis | Human, Mouse |
| ssGSEA | Single-sample GSEA | Human, Mouse |
|  | **---IMMUNE REPERTOIRE---** |  |
| TRUST4 | TCR and BCR Sequences Analysis | Human, Mouse |
|  | **---IMMUNE INFILTRATION---** |  |
| ImmuneDeconv | Cell Components Estimation | Human |
| mMCP | Immune Cell Estimation from mouse | Mouse |
|  | **---IMMUNE RESPONSE---** |  |
| MSIsensor2 | Microsatellite Instability (MSI) Detection | Human |
|  | **---MICROBIOME---** |  |
| Centrifuge | Bacterial Abundance Detection | Human, Mouse |
| PathSeq | Microbial sequences Detection | Human |
|  | **---NEO-ANTIGEN---** |  |
| arcasHLA | HLA Class I and Class II Genes Typing | Human |
|  | **---RIMA REPORT---** |  |
| report | RIMA HTML Report Using Multiqc | Human, Mouse |



## 1. Install

You can not run any jobs on the login node, even conda install is not allowed.


Download RIMA_pipeline folder to your own working directory.

### For human data:
```
git clone https://github.com/liulab-dfci/RIMA_pipeline.git
```

### For mouse data: 
```
## In the RIMA_pipeline folder, change to the RIMA-mouse branch

cd ./RIMA_pipeline
git checkout RIMA_Mouse
```

## 2. Link the reference folder

### For human reference (hg38)
```
## In the RIMA_pipeline folder
## This command will create a symbolic link to the human reference on Kraken

ln -s /data1-common/RIMA_references/ref_files
```

### For mouse reference (mm10)
```
## In the RIMA_pipeline folder
## This command will create a symbolic link to the mouse reference on Kraken

ln -s /data1-common/RIMA_references/mm10_ref ref_files
```


## 3. Activate the RIMA enviroment

```
export CONDA_ROOT=/liulab/linyang/rnaseq_env/miniconda3
export PATH=/liulab/linyang/rnaseq_env/miniconda3/bin:$PATH

source activate /liulab/linyang/rnaseq_env/miniconda3/envs/rnaseq
```

## 4. Prepare the 2 required execution files (metasheet.txt, config.yaml)

### 4.1 Example of metasheet

Ensure your metasheet contains **Two Required Columns** (SampleName, PatName) in comma-delimited format.
You can also add more phenotype information that you may want to compare e.g. columns for Responder, Age, Sex etc.

```
SampleName,PatName,Responder,Age,Gender
SRR8281231,P3,R,41,Male
SRR8281224,P13,NR,55,Female
SRR8281221,P20,NR,63,Male
SRR8281223,P21,NR,59,Male
......
```

### 4.2 Example of config.yaml

In the rnaseq_pipeline folder, we have prepared a config.yaml for you. 

First, ensure the data info matches data in the **Data information** section:

```
#########Fixed and user-defined parameters################
metasheet: metasheet.csv  # Meta info 
ref: ref.yaml             # Reference config 
assembly: hg38
cancer_type: GBM          #TCGA cancer type abbreviations
rseqc_ref: house_keeping  #Option: 'house_keeping' or 'false'. 
                          #By default, a subset of housekeeping genes is used by RSeQC to assess alignment quality.  
                          #This reduces the amount of time needed to run RSeQC.  
mate: [1,2]               #paired-end fastq format, we recommend naming paired-end reads with _1.fq.gz and _2.fq.gz


#########Cohort level analysis parameters################
design: Group             # Condition on which to do comparsion (as set up in metasheet.csv)
Treatment: R              # Treatment use in DESeq2, corresponding to positive log fold change
Control: NR               # Control use in DESeq2, corresponding to negative log fold change
batch: syn_batch          # Options: 'false' or a column name from the metasheet.csv.  
                          # If set to a column name in the metasheet.csv, the column name will be used for batch effect analysis (limma)
                          # It will also be used as a covariate for differential analysis (DESeq2) to account for batch effect.  

pre_treated: false        # Option: true or false. 
                          # If set to false, patients are treatment naive.  
                          # If set to true, patients have received some form of therapy prior to the current study.

```

Parameters within square brackets should be updated to match your analysis goals. In this tutorial, **[Responder]** is the phenotype of interest for comparison as specified in the **metasheet.txt**. All the comparison figures stored under **/images** folder

```
############################################################
#                     list samples                         #
############################################################

samples:
  SRR8281228:
    - /mnt/zhao_trial/Zhao2019_PD1_Glioblastoma_RNASeq/SRR8281228_1.fastq.gz
    - /mnt/zhao_trial/Zhao2019_PD1_Glioblastoma_RNASeq/SRR8281228_2.fastq.gz
  SRR8281231:
    - /mnt/zhao_trial/Zhao2019_PD1_Glioblastoma_RNASeq/SRR8281231_1.fastq.gz
    - /mnt/zhao_trial/Zhao2019_PD1_Glioblastoma_RNASeq/SRR8281231_2.fastq.gz
```

Finally, set the path of your data in the **list samples** section and set the number of runs for each sample (samples' name must be consistent with your metasheet.txt)

Currently, only **fastq files** are accepted as input (including fastq.gz).

## 5. Choose the tools you want to run

Use **execution.yaml** to control which tools to run in RIMA. Most downstream analysis require outputs from **DATA PROCESSING** module, so please run the **DATA PROCESSING** module first, then selecting which tools you want to use.

Example of execution.yaml:
```
## Note: Preprocess individual and cohort module necessary to get the alignment and quality results.
## Run the remaining modules only after these two modules.
preprocess_individual: true
preprocess_cohort: true

## Optional modules
## Note: The below modules are specialized modules, each dealing with specific targets.
## Make sure to run individual and cohort of each module to get all the results.

## Individual runs
immune_repertoire_individual: false
microbiome_individual: false


## Cohort runs
differential_expression_cohort: false
immune_infiltration_cohort: false
microbiome_cohort: false
....
```

## 6. Execution

### Step1: Check the pipeline with a dry run to ensure correct script and data usage.

```
snakemake -s rnaseq.snakefile -np
```
### Step2: submit the job.

After the dry-run success, please use sbatch to submit the job or run it on the working node:

### !!!Please not run RIMA on the interactive node!

```
#!/bin/bash
#SBATCH --job-name=RIMA
#SBATCH --mem=64G       # total memory need
#SBATCH -c 32 #number of cores

snakemake -s rnaseq.snakefile -k
```
**Note**: Argument -j that set the cores for parallelly run. (e.g. '-j 4' can run 4 jobs parallelly at the same time) 
**Note**: Argument -k that can skip the error independent run. (This argument can save lots of time for running data at the first time)

## 7. Output files

The output files from each tools are stored under `RIMA_kraken/analysis`

# How to Run RIMA on Kraken

## RIMA Workflow 
![](https://github.com/Lindky/RIMA_Kraken/blob/master/files/multiqc/Pipeline_Workflow_mqc.png)


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

**Notice**: When you run Mouse data, do not run the tools that aren't available for **Mouse Model**  (e.g. if you set ImmnueDeconv: true in the execution file, the results from ImmnueDeconv would be unreliable when you use Mouse data)


## 1.Install

You can not run any jobs on the login node, even conda install is not allowed.


Download RIMA_pipeline folder to your own working directory. Currently, RIMA_kraken supports hg38 and mm10 data
```
git clone https://github.com/Lindky/RIMA_Kraken.git
```

## 2.Activate the RIMA enviroment

```
export CONDA_ROOT=/liulab/linyang/rnaseq_env/miniconda3
export PATH=/liulab/linyang/rnaseq_env/miniconda3/bin:$PATH

source activate /liulab/linyang/rnaseq_env/miniconda3/envs/rnaseq
```

## 3. Prepare the 2 required execution files (metasheet.txt, config.yaml)

### 3.1 Example of metasheet

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

### 3.2 Example of config.yaml

In the rnaseq_pipeline folder, we have prepared a config.yaml for you. 

First, ensure the data info matches data in the **Data information** section:

```
---
############################################################
#                   Data information                       #
############################################################

ref: ref.yaml
assembly: hg38 #hg38 or mm10

cancer_type: GBM #short name of cancer type
metasheet: metasheet_latest.txt
designs: [Responder] #the column from metasheet which is used to do comparsion
############################################################
#                       level1                            #
############################################################

### star 
#Possible values are [ff-firststrand, ff-secondstrand, ff-unstranded, fr-firststrand, fr-secondstrand, fr-unstranded (default), transfrags]
library_type: 'fr-firststrand'
stranded: true


```
To see the details for each parameters: [parameters Interpretation](https://github.com/Lindky/RIMA_Kraken/blob/master/Parameters_description.md)

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
############################################################
#                       run settings                       #
############################################################
runs:
  run1:
    - SRR8281228
  run2:
    - SRR8281231

```

Finally, set the path of your data in the **list samples** section and set the number of runs for each sample (samples' name must be consistent with your metasheet.txt)

Currently, only **fastq files** are accepted as input (including fastq.gz).

## 3. Choose the tools you want to run

Use **execution.yaml** to control which tools to run in RIMA. Most downstream analysis require outputs from **DATA PROCESSING** module, so please run the **DATA PROCESSING** module first, then selecting which tools you want to use.

Example of execution.yaml:
```
##DATA PROCESSING
star: false
salmon: false
rseqc: false

##DIFFERENTIAL EXPRESSION
batch_removal: false
deseq2: true
gsea: false
ssgsea: false

....
```

## 4.Execution

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
#SBATCH --mem=32G       # total memory need
#SBATCH --time=96:00:00
#SBATCH -c 16 #number of cores


snakemake -s rnaseq.snakefile
```
**note**: Argument -j that set the cores for parallelly run. (e.g. '-j 4' can run 4 jobs parallelly at the same time) 

## 5.Output files

The raw output files from each tools are stored under `RIMA_kraken/analysis`

The processed output files (esed for visualization) are stored under `RIMA_kraken/files`

The output figures are stored under `RIMA_kraken/images`

**RIMA html report** will be stored under `RIMA_kraken/files`


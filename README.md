# How to Run RIMA on Kraken

## RIMA Workflow 
![](https://github.com/Lindky/RIMA_Kraken/blob/master/files/multiqc/Pipeline_Workflow_mqc.png)

### The full tutorial 
https://liulab-dfci.github.io/RIMA/

### The tutorial for running RIMA on Kraken
https://docs.google.com/document/d/1bdC-rRdIsmFszJ0PdsZROjFOS5BV-k6W1husbRSyRKQ/edit#

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




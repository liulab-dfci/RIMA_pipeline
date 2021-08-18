### We have prepared the RIMA reference folder with the current version we used. You can directly download from Iris server:

```
# We have nearly 69G reference files

wget http://cistrome.org/~lyang/ref.tar.gz 

```


# Steps to make your own custom reference files for RIMA pipeline

### Reference fasta

The human GDC hg38 fasta file is downloaded from [GDC website](https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files). 

### Gene annotation file (gtf)

The human gtf annotation file is downloaded from [GENCODE website](https://www.gencodegenes.org/human/). The current annotation file we used is V27.


### build STAR index

```bash
conda activate rna

## STAR Version: STAR_2.6.1d
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ./ref_files/v27_index --genomeFastaFiles GRCh38.d1.vd1.CIDC.fa --sjdbGTFfile gencode.v27.annotation.gtf
...
00:04:54 ..... started STAR run
00:04:54 ... starting to generate Genome files
00:05:57 ... starting to sort Suffix Array. This may take a long time...
00:06:11 ... sorting Suffix Array chunks and saving them to disk...
00:17:43 ... loading chunks from disk, packing SA...
00:19:31 ... finished generating suffix array
00:19:31 ... generating Suffix Array index
00:23:20 ... completed Suffix Array index
00:23:20 ..... processing annotations GTF
00:23:35 ..... inserting junctions into the genome indices
00:26:49 ... writing Genome to disk ...
00:27:06 ... writing Suffix Array to disk ...
00:28:53 ... writing SAindex to disk
00:29:05 ..... finished successfully
```
### RSeQC reference files

We download the human annotation bed file including the whole genome bed file, and house keeping bed file from RSeQC page from [sourcforge website](https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/).

```bash
./ref_files/refseqGenes.bed
./ref_files/housekeeping_refseqGenes.bed
```

### build salmon index

```bash
conda activate rna

## salmon Version: salmon 1.1.0
salmon index -t GRCh38.d1.vd1.CIDC.fa -i salmon_index

...
index ["salmon_index"] did not previously exist  . . . creating it
[jLog] [info] building index 
[jointLog] [info] [Step 1 of 4] : counting k-mers
[jointLog] [info] Replaced 164,553,847 non-ATCG nucleotides
[jointLog] [info] Clipped poly-A tails from 0 transcripts
[jointLog] [info] Building rank-select dictionary and saving to disk
[jointLog] [info] done
Elapsed time: 0.191866s
[jointLog] [info] Writing sequence data to file . . .
[jointLog] [info] done
Elapsed time: 1.91244s
[jointLog] [info] Building 64-bit suffix array (length of generalized text is 3,088,286,426)
[jointLog] [info] Building suffix array . . .
success
saving to disk . . . done
Elapsed time: 18.3072s
done
Elapsed time: 703.843s
```
### GMT file for gene set analysis

The GMT file is downloaded from [BROAD release page](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/6.1/). The current GMT file we used is "c2.cp.kegg.v6.1.symbols.gmt"


### STAR-Fusion genome resource lib

The genome resource lib is downloaded from [BROAD release page](https://www.gencodegenes.org/human/). The current lib we used is GRCh38_v22_CTAT_lib.

You can also prep it for use with STAR-fusion.
More details, read:

*  https://github.com/STAR-Fusion/STAR-Fusion/wiki/installing-star-fusion


### Centrifuge index

The human Centrifuge index is downloaded from [Centrifuge website](http://www.ccb.jhu.edu/software/centrifuge/). The current index we used is p_compressed+h+v that includes human genome, prokaryotic genomes, and viral genomes.

You can also build your own custom Centrifuge index. 
More details, read:

*  https://github.com/DaehwanKimLab/centrifuge

### TRUST4 reference files

TRUST4 reference files includes 1. TCR, BCR genomic sequence fasta file; 2. Reference database sequence containing annotation information.

```
hg38_bcrtcr.fa
human_IMGT+C.fa
```
These reference files can directlt be downloaded from [TRUST4 github](https://github.com/liulab-dfci/TRUST4).



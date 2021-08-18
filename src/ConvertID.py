#!/usr/bin/env python
"""Script to convert ensemble ID or transcript ID to gene symbol

"""
import pandas as pd
import sys
import os
from optparse import OptionParser


usage = "USAGE: %prog -f rsem.matrix"
optparser = OptionParser(usage=usage)
optparser.add_option("-m", "--matirx", help="file to process")
optparser.add_option("-r", "--reference", help="annotation file")
optparser.add_option("-p", "--prefix", help="output prefix")
optparser.add_option("-t", "--type", help="output type")
(options, args) = optparser.parse_args(sys.argv)

if not options.matirx:
  optparser.print_help()
  sys.exit(-1)


def Ensembl2name(df,colname,anno):
  ###convert ensembl id to gene symbol
  annot_dict = anno.set_index('ensembl_id').to_dict()['gene_name']
  df[colname] = df[colname].apply(lambda x: annot_dict[x])
  cols = df.shape[1]
  ###average the ensembl id with the same gene symbol
  df_redup = df.groupby('Gene_ID', as_index=False)[list(df.columns[1:cols])].mean() 
  return df_redup


def tx2name(df, colname, anno):
  ##convert transcript id to gene ensembl
  annot_dict = anno.set_index('transcript_id').to_dict()['ensembl_id']
  df["Gene_ID"] = df["Gene_ID"].apply(lambda x: annot_dict[x])
  cols = df.shape[1]
  ###sum transcripts of one gene ensembl
  df_sum = df.groupby('Gene_ID', as_index=False)[list(df.columns[1:cols])].sum()
  ###convert ensembl id to gene symbol
  df_redup = Ensembl2name(df_sum, "Gene_ID", anno)
  return df_redup



def name_entrezID(df, colname, anno):
  ##convert gene symbol to entrize id
  annot_dict = anno.set_index('gene_symbol').to_dict()['entrez_id']
  df["gene_id"] = df["gene_id"].map(annot_dict).fillna(df["gene_id"])
  cols = df.shape[1]
  transfered_redup = df.groupby('gene_id', as_index=False)[list(df.columns[1:(cols)])].mean()
  df_redup = transfered_redup[transfered_redup['gene_id'].isin(list(annot_dict.values()))]
  return df_redup

    

if __name__=='__main__':
  data = pd.read_csv(options.matirx)
  anno = pd.read_csv(options.reference, delimiter='\t',header=None,names=["transcript_id","ensembl_id","gene_name"])
  #####convert ID
  if options.type == "rsem":
    data_redup = Ensembl2name(data, "Gene_ID", anno)
  if options.type == "salmon":
    data_redup = tx2name(data, "Gene_ID", anno) 
  if options.type == "tide":
    anno = pd.read_csv(options.reference,sep = "\t" ,low_memory=False, names=["gene_symbol","entrez_id"], dtype = {'entrez_id': str})
    data_redup = name_entrezID(data, "gene_id", anno)
  data_redup.to_csv(options.prefix,index = False,sep = "\t")


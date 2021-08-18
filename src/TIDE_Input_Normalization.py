#!/usr/bin/env python
"""Script to normalize tumor sample expression by averaging normal controls
"""

import os
import pandas as pd
import numpy as np
import sys
from optparse import OptionParser


def main():
    usage = "USAGE: %prog -f rsem.matrix"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-i", "--input", help="expression file")
    optparser.add_option("-m", "--metadata", help="meta information for samples")
    optparser.add_option("-p", "--prefix", help="output prefix")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.input:
        optparser.print_help()
        sys.exit(-1)

    ###read in data
    expr = pd.read_csv(options.input, sep = ",", header = 0)
    metadata = pd.read_csv(options.metadata, sep = ',')


    ##take normal samples as control
    if len(set(metadata["Tissue"].tolist())) >1 : 
        ###average expression of normal controls
        t_all = list(map(''.join, itertools.product(*(sorted(set((c.upper(), c.lower()))) for c in "tumor"))))
        n_all = list(map(''.join, itertools.product(*(sorted(set((c.upper(), c.lower()))) for c in "normal"))))
        controls = metadata.loc[metadata["Tissue"] in n_all, "SampleName"]
        tumors = metadata.loc[metadata["Tissue"] in t_all, "SampleName"]
        expr["aver_control"] = expr[controls].mean(axis = 1)
        expr_normalize = expr[tumors].sub(expr["aver_control"], axis = 0)
    ##take the average of all samples as control
    else: 
        all_samples = metadata.loc[:, "SampleName"]
        expr_normalize = expr[all_samples].sub(expr[all_samples].mean(axis = 1), axis=0) #df.sub(df.mean(axis=1), axis=0)

    ###insert gene id in the first column
    expr_normalize.insert(0, column = "gene_id", value = expr["gene_id"])
    expr_normalize.to_csv(options.prefix+'_normalize_control.txt',index = False, sep = "\t")

if __name__ == '__main__':
	main()
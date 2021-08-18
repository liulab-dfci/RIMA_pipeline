#!/usr/bin/env python

#Script to make config file for fastq sample

import os
import sys
import pandas as pd

from optparse import OptionParser

def main():
    usage = "USAGE: %prog -m [meta] -o [output file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-m", "--meta", help="meta info")
    optparser.add_option("-o", "--out", help="output file")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.meta or not options.out:
        optparser.print_help()
        sys.exit(-1)

    out = open(options.out, "w")
    out.write("ref: ref.yaml\n")
    out.write("assembly: mm10\n")
    out.write("metasheet: metasheet.csv\n")
    out.write("designs: [Location]\n")
    out.write("\n")

    out.write("### star\n")
    out.write("library_type: 'fr-firststrand'\n")
    out.write("threads: 8\n")
    out.write("stranded: ture\n")

    out.write("\n")

    out.write("### rseqc\n")
    out.write("rseqc_ref: house_keeping\n")
    out.write("\n")

    out.write("### batch removal\n")
    out.write("batch_covariates: [no]\n")
    out.write("\n")

    out.write("### deseq2\n")
    out.write("batch: [no]\n")
    out.write("comparison: between\n")
    out.write("\n")

    out.write("### microbia\n")
    out.write("centrifuge: true\n")
    out.write("\n")

    out.write("samples:\n")
    meta = pd.read_csv(options.meta, sep=",")
    ID = meta['SampleName']
    fastq = meta['Fastq']

    for f in range(len(fastq.index)):
        tmp_ID = ID[f]
        tmp_fastq = fastq[f].split(';')
        out.write("  %s:\n" % ID[f])
        if len(tmp_fastq) == 2:
            out.write("    - %s\n" % tmp_fastq[0])
            out.write("    - %s\n" % tmp_fastq[1])
        else:
            out.write("    - %s\n" % tmp_fastq[0])
            
    out.close()

if __name__=='__main__':
    main()

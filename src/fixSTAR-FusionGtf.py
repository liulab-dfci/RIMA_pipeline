#!/usr/bin/python
"""Fix the STAR-Fusion gtf by finding missing transcript ids and setting
them to gene_id"""
import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog [options] file1 file2 ... fileN"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--file", help="GTF")

    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.file:
        optparser.print_help()
        sys.exit(-1)
    
    f = open(options.file)
    for l in f:
        if l.startswith("#"):
            continue
        tmp = l.strip().split("\t")
        #transcript attributes are the last column, key-value separated by ;
        attribs = (tmp[-1]).split(";")
        #for convenience we store it in a dict--
        #NOTE: there's a trailing semi-colon so we have to check for non-Null
        attrib_dic= dict([tuple(elm.split()) for elm in attribs if elm])
        
        #FINALLY we check for transcript_id
        if "transcript_id" not in attrib_dic:
            #insert into end of the list
            attrib_dic["transcript_id"] = attrib_dic["gene_id"]

        #compose the attribs list again
        a = "; ".join(["%s %s" % (k,v) for k,v in attrib_dic.items()])
        #substitute attribs back into the line
        tmp[-1] = a
        print("\t".join(tmp))
        
if __name__=='__main__':
    main()


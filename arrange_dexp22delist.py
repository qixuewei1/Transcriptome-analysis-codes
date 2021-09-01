#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
import argparse
import subprocess
import sys
from collections import defaultdict
import pandas as pd


def get_deg_dict(infile):
    diffsets = set()
    diffData = pd.read_csv(infile,header=0,sep="\t",index_col=0)
    upData = diffData.loc[(diffData.FDR<0.05) & (diffData.log2FC >1)]
    downData = diffData.loc[(diffData.FDR<0.05) & (diffData.log2FC <-1)]
    upset = set(upData.index)
    downset = set(downData.index)
    diffset = upset | downset
    diffsets = diffsets | diffset
    return diffsets

parser = argparse.ArgumentParser(description="huoqu can shu")
parser.add_argument('-i', '--exp',required=True, help=' file gff3')
parser.add_argument('-o', '--odr',required=True, help='ouput file fasta')
parser.add_argument('-f', '--ouf',required=True, help='ouput file name')
args = parser.parse_args()
outdir = os.path.abspath(args.odr)
print outdir
with open(args.exp)as exp,open(outdir+"/"+args.ouf,"w")as ogo:
    for line in exp:
        defile = line.strip("\n")
        namefile = defile.split("/")[-1]
        diffset = get_deg_dict(defile)
        with open(outdir+"/"+namefile+"_list","w")as dell:
            for eachgene in diffset:
                dell.write(eachgene+"\n")
        ogo.write(outdir+"/"+namefile+"_list\n")



# -*- coding: utf-8 -*-
import argparse
import pandas as pd
import numpy as np
import os,sys


if __name__ == "__main__" :
    parser = argparse.ArgumentParser(description="qu can")
    parser.add_argument( '--grouplist',required=True, help='diffience gene list')
    parser.add_argument( '--compair',required=True, help='diffience gene list')
    parser.add_argument( '--outdir',required=True, help='diffience gene list')
    args = parser.parse_args()
    outdir = os.path.abspath(args.outdir)
    file_list = []
    groupmatrix = pd.read_csv(args.grouplist,header=0,sep="\t",index_col=None,names=["sampleid","groupid"])
    with open(args.compair,"r")as compair:
        n = 0 
        for line in compair:
            n += 1
            if n !=1:
                group_list = line.strip("\n").split("\t")
                namef = "_".join(group_list)
                getmatrix = groupmatrix.loc[(groupmatrix["groupid"]==group_list[0]) | (groupmatrix["groupid"]==group_list[1])]
                getmatrix.to_csv(outdir+"/"+namef+"_group.list",header=True,sep="\t",index=None)
                file_list.append(outdir+"/"+namef+"_group.list")
    indexfile = outdir+"/group.list.list"
    with open(indexfile,"w")as indexf:
        outstr = "\n".join(file_list)
        indexf.write(outstr)




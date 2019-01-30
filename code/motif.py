#!/usr/bin/env python
#_*_coding:utf-8_*_


import re
import sys

def motif(fasta):
        l=len(fasta)
        motifs=[];
        b="X" * int(8)
        c=b+fasta+b
        n = 17
        li = [ c[i:i+n] for i in range(len(c)-n+1) ]	## Creating pattern of length 17 #######
        motifs.extend(li);
        motifs1=pd.DataFrame(motifs)
        return motifs1;

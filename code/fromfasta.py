#!/usr/bin/env python
#_*_coding:utf-8_*_


## This converts fasta to normal string sequence
import pandas as pd
import numpy as np

def fromfasta(file):
    df = pd.read_csv(file,header=None, sep='\n',error_bad_lines=False);
    newdf = [];
    for i in range(len(df)):
        if(df.iloc[i][0][0]!='>'):
            newdf.append(df.iloc[i][0]);
    newdf=pd.DataFrame(newdf)
    return newdf

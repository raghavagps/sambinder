#!/usr/bin/env python
#_*_coding:utf-8_*_

#This code takes Fasta Sequence and changes them to normal dataframe with the sequence names

import pandas as pd

def fromfasta(file):
    df = pd.read_csv(file,header=None, sep='\n',error_bad_lines=False);
    newdf = [];
    temp='';
    l=len(df)
    sequencenames=[];
    for i in range(len(df)):
        
        if(df.iloc[i][0][0]!='>' or i==l-1):
            
            temp=temp+str(df.iloc[i][0]);
            if(i==l-1):
                newdf.append(temp);
                sequencenames.append(df.iloc[i][0]);
        
        elif(df.iloc[i][0][0]=='>' and i!=l-1):
            if(i!=0):
                newdf.append(temp);
                temp=''
            sequencenames.append(df.iloc[i][0]);


    newdf=pd.DataFrame(newdf)
    return newdf,sequencenames[:-1]

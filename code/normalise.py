#!/usr/bin/env python
#_*_coding:utf-8_*_


import sys
import math
def normalise(df):
    
    #print('In Normalise(df), printing df \n',df.iloc[:][0])

    minimum = df.min();
    
    maximum = df.max();
    
    rows = df.shape[0];
    norm=[0]*rows;
    #print(maximum,minimum)

    for i in range(rows):
        norm[i] = ((float(df.iloc[i][:])-minimum)/(maximum-minimum)*9.0)
        norm[i]=round(norm[i],0) 
        
    #norm1=pd.DataFrame(norm);
    #print(norm1)
    return norm;

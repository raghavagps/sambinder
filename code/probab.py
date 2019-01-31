#!/usr/bin/env python
#_*_coding:utf-8_*_

import pandas as pd

def probab(df):
    
    rows = df.shape[0];
    prob_val=[0.]*rows;
    for i in range(rows):
        prob_val[i] = ((float(df.iloc[i][:])*float(1.45))-float(0.20))
    prob_val=pd.DataFrame(prob_val);
      
    return prob_val;

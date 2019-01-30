
# coding: utf-8

# In[ ]:


## This converts fasta to normal string sequence
import pandas as pd
import numpy as np
import os
import re
import sys
import math
import warnings
warnings.filterwarnings("ignore")
from sklearn.svm import SVC
from sklearn.datasets import load_svmlight_file
from sklearn.externals import joblib
from sklearn.ensemble import RandomForestClassifier
import csv
import itertools

currdir=os.getcwd();

def fromfasta(file):
    df = pd.read_csv(file,header=None, sep='\n',error_bad_lines=False);
    newdf = [];
    temp='';
    for i in range(len(df)):
        if(df.iloc[i][0][0]!='>'):
            #temp=temp+df.iloc[i][0];
            newdf.append(df.iloc[i][0]);
        '''else:
                                    newdf.append(temp);
                                    temp='';
                            '''
    newdf=pd.DataFrame(newdf)
    print(newdf)
    return newdf


# In[ ]:


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


# In[ ]:


assign = {'A':[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
          'C':[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
          'D':[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
          'E':[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
          'F':[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
          'G':[0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
          'H':[0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
          'I':[0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
          'K':[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
          'L':[0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
          'M':[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
          'N':[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
          'P':[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
          'Q':[0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
          'R':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
          'S':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
          'T':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
          'V':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
          'W':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
          'Y':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
          'X':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}

def motif2bin(seq1):
    seq=pd.DataFrame(seq1)
    l=len(seq);
    #print(seq);
    
    #print(type(seq));
    binary=[]; # To store binary profile of all the motifs
    for i in range(0,l):
    
        temp=[]; # To store binary profile of each motif
        
        for j in range(0,len(seq.iloc[i][0])):
            
            #print(len(seq.iloc[i][0]))
                
            temp.extend(assign[seq.iloc[i][0][j]])
            #print('HELLO',j,i,seq.iloc[i][0][j],assign[seq.iloc[i][0][j]]);
        
        binary.append(temp);
        #binary1=pd.DataFrame(binary);
    
    return binary;


# In[ ]:

def probab(df):
    print('Entering probab()')
    rows = df.shape[0];
    prob_val=[0.]*rows;
    for i in range(rows):
        prob_val[i] = ((float(df.iloc[i][:])*float(1.45))-float(0.20))
    prob_val=pd.DataFrame(prob_val);
      
    return prob_val;

def normalise(df):
    
    print('Entering normalise()')
    #print('In Normalise(df), printing df \n',df.iloc[:][0])

    minimum = df.min();
    
    maximum = df.max();
    
    rows = df.shape[0];
    norm=[0.]*rows;
    print(maximum,minimum)

    for i in range(rows):
        print(float(df.iloc[i][:]));
        norm[i] = ((float(df.iloc[i][:])-minimum)/(maximum-minimum)*9.0)
        #norm[i]=round(norm[i],0) 
        
    #norm1=pd.DataFrame(norm);
    print(norm)
    return norm;


# In[ ]:


def predict(testdata, method,threshold):
    
    print('Entering predict()')
    df = pd.DataFrame();
    y_p_score1=pd.DataFrame();
    data_test = pd.DataFrame(testdata);
    X_test = data_test.iloc[:][0:356];
    a=[];
    k=0;
    
    if(method==1):
        path=currdir+'/Trained Models/svc_binary_model';
        clf = joblib.load(path)
        y_p_score1=clf.decision_function(X_test)
        #print(y_p_score1);
        
    elif(method==2):
        path=currdir+'/Trained Models/rf_binary_model';
        clf = joblib.load(path)
        y_p_score=clf.predict_proba(X_test)
        y_p_score1 = [i[1] for i in y_p_score] 
        
    elif(method==3):
        path=currdir+'/Trained Models/ann_binary_model';
        clf = joblib.load(path)
        y_p_score=clf.predict_proba(X_test)
        y_p_score1 = [i[1] for i in y_p_score] 
        
    elif(method==4):
        path=currdir+'/Trained Models/svc_pssm_model';
        clf = joblib.load(path)
        y_p_score1=clf.decision_function(X_test)
    
    elif(method==5):
        path=currdir+'/Trained Models/rf_pssm_model';
        clf = joblib.load(path)
        y_p_score=clf.predict_proba(X_test)
        y_p_score1 = [i[1] for i in y_p_score] 
    print(path)    
    interact =[];

    #y_p_score_temp = probab(pd.DataFrame(y_p_score1));

    y_p_score2=normalise(pd.DataFrame(y_p_score1));
    
    y_p_score = pd.DataFrame(y_p_score2);
    y_p_score=y_p_score.astype(int);
    for i in range(len(y_p_score)):
        if(int(y_p_score.iloc[i][0])<threshold):
            interact.extend('-');
        else:
            interact.extend('+');  
    print(y_p_score,interact)     
    return y_p_score,interact;c


# In[ ]:


def sambinder(mode):
    
    met = int(input('Choose desired prediction method: \n 1) Binary SVC \n 2) Binary Random Forest \n 3) Binary MLP \n 4) PSSM SVC \n 5) PSSM Random Forest \n\n \t'))
    
    if met<=0 or met>=6:
        print('\tNot a valid entry. Kindly choose a value between 1-5');
        sambinder();
    
    thresh = int(input('\n Enter the threshold probability score (between 0 to 9) \n\n\t'));
    
    if thresh<0 or thresh >9:
        print('\tNot a valid threshold');
        sambinder();
    path = input('\n Enter the path to your FASTA sequence\n\n\t')
        
    if(type(path)!=str):
        print("Enter path to your FASTA sequence");
        
        
        
    nofasta = fromfasta(path);
    
    
    nofasta = nofasta.values.tolist();
    
    nofasta = list(itertools.chain.from_iterable(nofasta))
    #print('LIST\n',nofasta)
    output=[];
    
    for i in range(len(nofasta)):
            
        #print("#####", i,"####")
        
        if(met==1 or met==2 or met==3):

            
            motifs = motif(nofasta[i]);
            
            
            binary = motif2bin(motifs);
            
            score,interact = predict(binary,met,thresh);
            
            #print(type(motifs),type(binary),type(score),type(interact))
            
            
            score = score.values.tolist();
            
            output.append([">seq_"+str(i+1)]);
            
            #print("Appending 1st line:\n",nofasta[i][0])
            nofasta[i] = list(itertools.chain.from_iterable(nofasta[i]))
            output.append(nofasta[i]);
            
            #print("Appending 2nd line:\n",score)
            score = list(itertools.chain.from_iterable(score))
            output.append(score);
            
            #print("Appending 3rd line:\n",interact)
            
            output.append(interact);
            
            #print("Printing Output\n",output)
            
            #print(output)
    #output1=pd.DataFrame(output);
    
    #csvout=output1.to_csv();
    file = open('SamBinderOut.csv','w')
    with file:
        writer = csv.writer(file);
        writer.writerows(output);
    return output;    


# In[ ]:


sambinder('abc')


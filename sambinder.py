################################# This code is developed to predict SAM interacting residue in the protein chain ######
################################# Developed by G.P.S. Raghava, IIIT Delhi, New Delhi, India ###########################
################################################### Date 30-01-2019 ###################################################
 

# coding: utf-8

# In[1]:


## This converts fasta to normal string sequence
import pandas as pd
import numpy as np
import os
import re
import sys
import math
from sklearn.svm import SVC
from sklearn.datasets import load_svmlight_file
from sklearn.externals import joblib
from sklearn.ensemble import RandomForestClassifier
import csv
import itertools
import getopt
import warnings
warnings.filterwarnings("ignore")
currdir=os.getcwd();

###################################### Function for Reading Inputfile in FASTA format ##############################
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

# In[2]:
######################################### Function for Generating patterns of fixed length 17 ####################### 
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

# In[3]:

######################################### Function for Assigning Binary values and Generating Binary Profile ####################### 

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
    
    binary=[]; # To store binary profile of all the motifs
    for i in range(0,l):
    
        temp=[]; # To store binary profile of each motif
        
        for j in range(0,len(seq.iloc[i][0])):
            temp.extend(assign[seq.iloc[i][0][j]])
        
        binary.append(temp);
    return binary;

# In[4]:

########################## Function for fitting Linear regression equation into the output machine learning score ################

def probab(df):
    #print('Entering probab()')
    rows = df.shape[0];
    prob_val=[0.]*rows;
    for i in range(rows):
        prob_val[i] = ((float(df.iloc[i][:])*float(1.45))-float(0.20))
    prob_val=pd.DataFrame(prob_val);
      
    return prob_val;

# In[5]:

######################################### Function for normalizing score into the propensity value in between 0-9 ####################### 

def normalise(df):
    
    minimum = df.min();
    maximum = df.max();
    
    rows = df.shape[0];
    norm=[0.]*rows;

    for i in range(rows):
        #print(float(df.iloc[i][:]));
        norm[i] = ((float(df.iloc[i][:])-minimum)/(maximum-minimum)*9.0)
        norm[i]=round(norm[i],0)
        
    return norm;


# In[ ]:

######################################### Function for Normalizing PSSM matrix  ####################### 

def pssm_n1(file):
    
    filename,file_ext=os.path.splitext(file)
    df=pd.read_csv(file,header=None,error_bad_lines=False)
    df1 = df.iloc[:,0:22]
    def pssm1(x):
        if type(x) is str:
            return x
        elif x:
            return (1/(1+(2.7182)**(-x)))
        else:
            return
    df2 = df1.applymap(pssm1)
    
    df3 = df2.round(2)
    filenamee = filename+"_pssm_n1.csv"
    df3.to_csv(filenamee, encoding='utf-8', index=False, header=False)
    return df3,filenamee

# In[ ]:

######################################### Function for Generating PSSM profiles  ####################### 

import sys, getopt
import numpy
import itertools

def pssm2pat(inputfile):
    #inputfile = ''
    filename,file_ext=os.path.splitext(inputfile)
    outputfile = filename+'_Pattern.csv';
    window= 17
    with open(inputfile,'r') as f:
        g = list(f)
    
    orig_stdout = sys.stdout
    n = open(outputfile,'w')
    sys.stdout = n
    no_insert_seq = (window-1)/2
    x1 = "X, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0\n"
    i = 0
    j = 0
    i1 = 0
    j1 = 0
    l = 0
    b = []
    
    while j < no_insert_seq:
        g.insert(0,x1)
        g.insert(len(g),x1)
        j += 1
    
    while i < len(g) :
        g[i] = g[i].replace('\n','')
        g[i] = g[i].replace("'",'')
        g[i] = g[i].split(",")
        i += 1
    hh = numpy.zeros(shape=(len(g),len(g[0])-1),dtype=float)
    while i1 < len(g) :
        while j1 < len(g[0]) :
            if j1 > 0 :
                if(type(g[i1][j1])!=str):
                    hh[i1][j1-1] =float(g[i1][j1]);
                else:
                    hh[i1][j1-1] = 0.0;
            j1 += 1
        j1 = 0
        i1 += 1
    
    while l < len(hh) :
        if (l + window) < len(g)+1 :
            b = hh[l:l+window]
            print(str(list(itertools.chain(*b))).replace('[','').replace(']','').replace("'",""))
        l += 1
    n.truncate();


# In[ ]:

######################################### Function for perorming Machine Learning Prediction ####################### 

def predict(testdata, method,threshold):
    
    df = pd.DataFrame();
    y_p_score1=pd.DataFrame();
    data_test = pd.DataFrame(testdata);
    X_test = data_test.iloc[:][0:356];
    a=[];
    k=0;
    
    if(method==1):
        path=currdir+'/sam_models/svc_binary_model';
        clf = joblib.load(path)
        y_p_score1=clf.decision_function(X_test)
        
    elif(method==2):
        path=currdir+'/sam_models/rf_binary_model';
        clf = joblib.load(path)
        y_p_score=clf.predict_proba(X_test)
        y_p_score1 = [i[1] for i in y_p_score] 
        
    elif(method==3):
        path=currdir+'/sam_models/ann_binary_model';
        clf = joblib.load(path)
        y_p_score=clf.predict_proba(X_test)
        y_p_score1 = [i[1] for i in y_p_score]
 
    elif(method==4):
        path=currdir+'/sam_models/knn_binary_model';
        clf = joblib.load(path)
        y_p_score=clf.predict_proba(X_test)
        y_p_score1 = [i[1] for i in y_p_score]
    
    elif(method==5):
        path=currdir+'/sam_models/svc_pssm_model';
        clf = joblib.load(path)
        y_p_score1=clf.decision_function(X_test)
    
    elif(method==6):
        path=currdir+'/sam_models/rf_pssm_model';
        clf = joblib.load(path)
        y_p_score=clf.predict_proba(X_test)
        y_p_score1 = [i[1] for i in y_p_score] 
    interact ='';
    new_score='';
    
    y_p_score2=normalise(pd.DataFrame(y_p_score1));    
    y_p_score = pd.DataFrame(y_p_score2);
    y_p_score=y_p_score.astype(int);
    for i in range(len(y_p_score)):
        if(int(y_p_score.iloc[i][0])<threshold):
            interact=interact+'-';
        else:
            interact=interact+'+'; 
        new_score=new_score+str(y_p_score.iloc[i][0])
    return new_score,interact;

# In[8]:

######################################### SAMbinder Validation and Processing Function ####################### 

def sambinder(path,outputfile,thresh,met):
    
    inputfile, input_ext = os.path.splitext(path);
    if met<=0 or met>=6:
        print('\tNot a valid entry. Kindly choose a value between 1-5');
        sambinder();
    
    if thresh<0 or thresh >9:
        print('\tNot a valid threshold');
        sambinder();
    
    nofasta1,sequencenames = fromfasta(path);
    nofasta = nofasta1.values.tolist();
    nofasta = list(itertools.chain.from_iterable(nofasta))
    dfTemp=pd.DataFrame();
    output=[];
    for i in range(len(nofasta)):
        
        if(met==1 or met==2 or met==3 or met==4):
            motifs = motif(nofasta[i]);
            binary = motif2bin(motifs);
            score,interact = predict(binary,met,thresh);

            output.append([sequencenames[i]]);

            output.append([nofasta[i]]);
            l = len(nofasta[i]);
            
            output.append([score]);
            output.append([interact]);
    
        elif(met==5 or met==6):
            
            ans = nofasta[i]
            dfTemp=pd.DataFrame([nofasta1.iloc[i][:]]);
            dfTemp.to_csv('temp.fa',encoding='utf-8', index=False, header=False)
            
            codeDir = currdir+'/code/seq2pssm_imp.py';
            
            cmd = "python "+ codeDir+" -i temp.fa"+ " -o temp -d swissprot";
            
            os.system(cmd);
            df,filenamee = pssm_n1('PSSMProfile.txt');
            pssm2pat(filenamee);
            filenamee=filenamee[:-4];
            testdata=pd.read_csv(filenamee+'_Pattern.csv',header=None,error_bad_lines=False)
            message = 'rm temp* PSSMProfile_pssm_n1_Pattern.csv PSSMProfile_pssm_n1.csv PSSMProfile.txt';
            os.system(message);
            score,interact = predict(testdata,met,thresh);

            output.append([sequencenames[i]]);

            output.append([nofasta[i]]);
            l = len(nofasta[i]);
            
            output.append([score]);
            output.append([interact]);

    filename=outputfile+'.csv';
    file = open(filename,'w')
    with file:
        writer = csv.writer(file);
        writer.writerows(output);
    return output;

# In[9]:

######################################### SAMbinder Main Function ####################### 

def main(argv):
    
    inputfile = ''
    outputfile = ''
    
    if len(argv[1:]) == 0:
    
        print ("\nProgram to predict interacting and non-interacting sites in given peptide sequences\n")
        print ("Usage:\tpython3 sambinder.py -i <input file> -o <output file> -t <threshold> -m <method> \n")
      
        print ('-i\tInput file having sequences in FASTA format \n-o\toutputFile generated by SAMbinder\n-t\tThreshold for classifying a probability score as interacting and non-interacting \n-m\tMethod to be used\n \t\t 1) Binary SVC \n \t\t 2) Binary Random Forest \n\t\t 3) Binary MLP \n\t\t 4) Binary KNN \n\t\t 5) PSSM SVC \n\t\t 6) PSSM Random Forest \n\n')
        
        print('Example:\t python3 sambinder.py -i sample.fa -o sample -t 4 -m 2')
        
        sys.exit()

    try:
           # opts is a list of returning key-value pairs, args is the options left after striped
           # the short options 'hi:o:', if an option requires an input, it should be followed by a ":"
           # the long options 'ifile=' is an option that requires an input, followed by a "="
        opts, args = getopt.getopt(argv,"i:o:t:m:",["ifile=","ofile=","threshold=","method="])
    except getopt.GetoptError:
        print ('sambinder.py -i <input_file> -o <output_file> -t <threshold> -m <method>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('sambinder.py -i <input_file> -o <output_file> -t <threshold> -m <method>')
            sys.exit()
        elif opt in ("-i", "ifile"):
            inputfile = arg
        elif opt in ("-o", "ofile"):
            outputfile = arg
        elif opt in ("-t", "thresh"):
            thresh = int(sys.argv[6])
        elif opt in ("-m", "method"):
            method = int(arg)
    filename, file_extension = os.path.splitext(inputfile)

    sambinder(inputfile,outputfile,thresh,method);
main(sys.argv[1:])


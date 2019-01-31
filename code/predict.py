#!/usr/bin/env python
#_*_coding:utf-8_*_
from sklearn.svm import SVC
from sklearn.datasets import load_svmlight_file
from sklearn.externals import joblib
from sklearn.ensemble import RandomForestClassifier
import re
import numpy as np
import pandas as pd
import sys
import csv
import itertools

def predict(testdata, method,threshold):
    
    #print('Test data',testdata,'\n\n')
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
        #print(y_p_score1);
        
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

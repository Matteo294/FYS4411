# Common imports
import os, fnmatch
import csv
import re
import multiprocessing as mp
import pandas as pd
from pandas import DataFrame
import matplotlib
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
#from typing_extensions import Required
import numpy as np
import sys
import argparse
import time #to implement the running time 
from numpy import log2, zeros, mean, var, sum, loadtxt, arange, array, cumsum, dot, transpose, diagonal, sqrt
from numpy.linalg import inv



#instruction how to run the program 
##################################
#SELECTOR FROM COMMAND LINE
##################################
ap = argparse.ArgumentParser(description="The program performs the blocking analysis on a given bunch of data. User can select between different operating mode")
ap.add_argument("-p", "--Parallel_simulation", choices=['0', '1'], required=True, help=" 0--> standard_file_analysis; 1--> parallel_file_analysis")
ap.add_argument("-c", "--Program_selector", choices=['0','1','2','3'], required=True, help="Chose Program:  \
    \n 0--> Simple case one file analysis;\
    \n 1--> multiple file varying alpha; \
    \n 2--> multiple file varying dt;\
    \n 3--> multiple file varying N")
args = vars(ap.parse_args())

if args['Parallel_simulation']=='1': 
    dir = './Data/parallel'
elif args['Parallel_simulation']=='0':
    dir = "./Data/standard"

selector =int(args['Program_selector'])

#args = vars(ap.parse_args())
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# Choose the repository 
DATA_ID = "./"
def data_path(dat_id):
    return os.path.join(DATA_ID, dat_id)


#########################################
#definition of functions
#########################################
#Blocking alghoritm
def block(x):
    n = len(x)
    d = int(log2(n))
    s, gamma = zeros(d), zeros(d)
    mu = np.mean(x)
    var_to_compare = np.var(x)
    std_k0 = np.sqrt(var_to_compare/n)

    # estimate the auto-covariance and variances 
    # for each blocking transformation
    for i in arange(0,d):
        n = len(x)
        # estimate autocovariance of x
        gamma[i] = (n)**(-1)*sum( (x[0:(n-1)]-mu)*(x[1:n]-mu) )
        # estimate variance of x
        s[i] = np.var(x)
        # perform blocking transformation
        x = 0.5*(x[0::2] + x[1::2])
   
    # generate the test observator M_k from the theorem
    M = (cumsum( ((gamma/s)**2*2**arange(1,d+1)[::-1])[::-1] )  )[::-1]

    # we need a list of magic numbers
    q =array([6.634897,9.210340, 11.344867, 13.276704, 15.086272, 16.811894, 18.475307, 20.090235, 21.665994, 23.209251, 24.724970, 26.216967, 27.688250, 29.141238, 30.577914, 31.999927, 33.408664, 34.805306, 36.190869, 37.566235, 38.932173, 40.289360, 41.638398, 42.979820, 44.314105, 45.641683, 46.962942, 48.278236, 49.587884, 50.892181])

    # use magic to determine when we should have stopped blocking
    for k in arange(0,d):
        if(M[k] < q[k]):
            break
    if (k >= d-1):
        print("Warning: Use more data")
    
    for i in arange(0,d):
        s[i] = s[i]/2**(d-i)
    
    return mu, s[0:d], k, var_to_compare, std_k0


#simple run
def simplerun(dir):
    print ("\n=========================================\n")
    rep = "/singlerun/"
    filelist = sorted_alphanumeric(os.listdir(dir+rep))
    x=[]
    for f in filelist:
        infile = data_path(dir+rep+f)
        x = np.concatenate((x,np.genfromtxt(infile)))

    
    (mean, var, k, var_to_compare, std_k0) = block(x) 
    
    std = sqrt(var)
    data ={'Mean':[mean], 'STDev':[std[k]], 'Var_to_compare':[var_to_compare], 'std_k0':[std_k0]}
    print(data)
    print ("\n=========================================\n")


#Varying alpha 
def var_alpha(dir):
    rep="/varying_alpha/" 
    filelist = sorted_alphanumeric(fnmatch.filter(os.listdir(dir+rep), "*.dat"))
    alpha_val = np.genfromtxt('./Data/parallel/varying_alpha/varying_alpha.csv',skip_header=True)

    with open(os.path.join(dir+rep, "post_analysis_alpha"+'.csv'), "w", newline='') as fcsv:
        writer = csv.writer(fcsv,delimiter =',')
        writer.writerow(["alpha,energy","std"])
        for i in range(0,int(len(fnmatch.filter(os.listdir(dir+rep), '*core0*.dat')))): # the variable i cycles on the different alpha
            x=[]
            for f in filelist:   
                if fnmatch.fnmatch(f, '*alpha'+str(i)+'.dat'): # here we collect all the data from file that represent the same alpha 
                    infile = data_path(dir+rep+f)
                    x = np.concatenate((x,np.genfromtxt(infile)))
            
            (mean, var ,k, var_to_compare, std_k0) = block(x) #blocking analysis
            std = sqrt(var)
            data ={'Alpha':[alpha_val[i]], 'Mean':[mean], 'STDev':[std[k]], 'Var_to_compare':[var_to_compare], 'std_k0':[std_k0]}
            writer = csv.writer(fcsv,delimiter =',')
            writer.writerow([alpha_val[i],mean,std[k]])
            print(data)

    fcsv.close()

def var_dt(dir):
    rep= "/varying_dt/"
    filelist=sorted_alphanumeric(os.listdir(dir+rep))
    dt_val = np.genfromtxt('./Data/parallel/varying_dt/varying_dt.csv',skip_header=True)
    
    with open(os.path.join(dir+rep, "post_analysis_dt"+'.csv'), "w", newline='') as fcsv:
        writer = csv.writer(fcsv,delimiter =',')
        writer.writerow(["dt,energy","STD"])
        for i in range(0,int(len(fnmatch.filter(os.listdir(dir+rep), '*core0*.dat')))):
            x=[]
            for f in filelist:
                if fnmatch.fnmatch(f, '*dt'+str(i)+'.dat'):  
                    infile = data_path(dir+rep+f)
                    x = np.concatenate((x,np.genfromtxt(infile)))
            
            (mean, var ,k, var_to_compare, std_k0) = block(x) #blocking analysis
            std = sqrt(var)
            data ={'dt':[dt_val[i]],'Mean':[mean], 'STDev':[std[k]], 'Var_to_compare':[var_to_compare], 'std_k0':[std_k0]}
            writer = csv.writer(fcsv,delimiter =',')
            writer.writerow([dt_val[i], mean,std[k]])
            print(data)


def var_N(dir):
    rep= "/varying_N/"
    filelist =sorted_alphanumeric(os.listdir(dir+rep))
    N_val = np.genfromtxt('./Data/parallel/varying_N/varying_N.csv',skip_header=True)

    with open(os.path.join(dir+rep, "post_analysis_N"+'.csv'), "w", newline='') as fcsv:
        writer = csv.writer(fcsv,delimiter =',')
        writer.writerow(["N,energy","STD"])
        for i in range(0,int(len(fnmatch.filter(os.listdir(dir+rep), '*core0*.dat')))):
            x=[]
            for f in filelist:
                if fnmatch.fnmatch(f, '*N'+str(i)+'.dat'):  
                    infile = data_path(dir+rep+f)
                    x = np.concatenate((x,np.genfromtxt(infile)))
            (mean, var ,k, var_to_compare, std_k0) = block(x) #blocking analysis
            std = sqrt(var)
            data ={'N':[N_val[i]],'Mean':[mean], 'STDev':[std[k]], 'Var_to_compare':[var_to_compare], 'std_k0':[std_k0]}
            writer = csv.writer(fcsv,delimiter =',')
            writer.writerow([N_val[i], mean,std[k]])
            print(data)


#sort alphanurecaly the files
def sorted_alphanumeric(data):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(data, key=alphanum_key)


#################################################################
#MAIN
#################################################################
#Case 
if selector==0: #Simple run
    simplerun(dir)
if selector==1: #Varyng alpha
    print("Attention! The scripts analyzes all the scripts contained in the varying_alpha folder! If you just reduced the number of alphas you should manually delete files.")
    var_alpha(dir)
if selector==2: #Varyng dt
    print("Attention! The scripts analyzes all the scripts contained in the varying_dt folder! If you just reduced the number of alphas you should manually delete files.")
    var_dt(dir)
if selector==3: #Varyng N
    print("Attention! The scripts analyzes all the scripts contained in the varying_N folder! If you just reduced the number of alphas you should manually delete files.")
    var_N(dir)


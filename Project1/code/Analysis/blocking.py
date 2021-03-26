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
ap = argparse.ArgumentParser(description="The program performs the blocking analysis on a given bunch of data. User can select between different operating mode \
    and different bunches of data will be analyzed as a consequence. For each bunch of analyzed data the program prints on terminal \
    the mean value, the error estimated through the blocking procedure and the std/sqrt(N_MC_steps) evaluated on the whole set of data.")
ap.add_argument("-p", "--Parallel_simulation", choices=['0', '1'], required=True, help=" 0--> Analyzes data from ./standard; \
     1--> Analyzes data from ./parallel")
ap.add_argument("-c", "--Program_selector", choices=['0','1','2','3'], required=True, help="After choosing between ./parallel and ./standard, chose Program:  \
    0--> Analyzes data from prev_selection/singlerun/;\
    1--> Analyzes data from prev_selection/varying_alpha/; \
    2--> Analyzes data from prev_selection/varying_dt/;\
    3--> Analyzes data from prev_selection/varying_dt/ \
    Case 1,2,3 produce a .csv file with the resume of the analysis.")
args = vars(ap.parse_args())

if args['Parallel_simulation']=='1': 
    dir = './Data/parallel'
    seach_word_alpha = 'core0'
    seach_word_dt = 'core0'
    seach_word_N = 'core0'
elif args['Parallel_simulation']=='0':
    dir = "./Data/standard"
    seach_word_alpha = 'stepalpha'
    seach_word_dt = 'stepdt'
    seach_word_N = 'stepN'

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
    
    return mu, s[0:d], k


#simple run
def simplerun(dir):
    print ("\n=========================================\n")
    rep = "/singlerun/"
    filelist = sorted_alphanumeric(os.listdir(dir+rep))
    x=[]
    for f in filelist:
        infile = data_path(dir+rep+f)
        x = np.concatenate((x,np.genfromtxt(infile)))

    
    (mean, s, k) = block(x) 
    
    std = sqrt(s)
    data ={'Mean':[mean], 'sigma_block':[std[k]], 'full_std':std[0]}
    print(data)
    print ("\n=========================================\n")


#Varying alpha 
def var_alpha(dir):
    rep="/varying_alpha/" 
    filelist = sorted_alphanumeric(fnmatch.filter(os.listdir(dir+rep), "*.dat"))
    alpha_val = np.genfromtxt(dir+'/varying_alpha/varying_alpha.csv',skip_header=True)

    with open(os.path.join(dir+rep, "post_analysis_alpha"+'.csv'), "w", newline='') as fcsv:
        writer = csv.writer(fcsv,delimiter =',')
        writer.writerow(["alpha,energy,std"])
        for i in range(0, int(len(fnmatch.filter(os.listdir(dir+rep), '*'+ seach_word_alpha +'*.dat')))): # the variable i cycles on the different alpha
            x=[]
            for f in filelist:   
                if fnmatch.fnmatch(f, '*alpha'+str(i)+'.dat'): # here we collect all the data from file that represent the same alpha 
                    infile = data_path(dir+rep+f)
                    x = np.concatenate((x,np.genfromtxt(infile)))
        
            (mean, s ,k) = block(x) #blocking analysis
            std = sqrt(s)
            data ={'Alpha':[alpha_val[i]], 'Mean':[mean], 'signa_block':[std[k]], 'full_std':[std[0]]}
            writer = csv.writer(fcsv,delimiter =',')
            writer.writerow([alpha_val[i],mean,std[k]])
            print(data)

    fcsv.close()



def var_dt(dir):
    rep= "/varying_dt/"
    filelist=sorted_alphanumeric(os.listdir(dir+rep))
    dt_val = np.genfromtxt(dir+'/varying_dt/varying_dt.csv',skip_header=True)
    
    with open(os.path.join(dir+rep, "post_analysis_dt"+'.csv'), "w", newline='') as fcsv:
        writer = csv.writer(fcsv,delimiter =',')
        writer.writerow(["dt,energy,std"])
        for i in range(0,int(len(fnmatch.filter(os.listdir(dir+rep), '*'+seach_word_dt+'*.dat')))):
            x=[]
            for f in filelist:
                if fnmatch.fnmatch(f, '*dt'+str(i)+'.dat'):  
                    infile = data_path(dir+rep+f)
                    x = np.concatenate((x,np.genfromtxt(infile)))
            
            (mean, s ,k) = block(x) #blocking analysis
            std = sqrt(s)
            data ={'dt':[dt_val[i]],'Mean':[mean], 'sigma_block':[std[k]], 'full_std':[std[0]]}
            writer = csv.writer(fcsv,delimiter =',')
            writer.writerow([dt_val[i], mean,std[k]])
            print(data)
    
    fcsv.close()


def var_N(dir):
    rep= "/varying_N/"
    filelist =sorted_alphanumeric(os.listdir(dir+rep))
    N_val = np.genfromtxt(dir+'/varying_N/varying_N.csv',skip_header=True)

    with open(os.path.join(dir+rep, "post_analysis_N"+'.csv'), "w", newline='') as fcsv:
        writer = csv.writer(fcsv,delimiter =',')
        writer.writerow(["N,energy,std"])
        for i in range(0,int(len(fnmatch.filter(os.listdir(dir+rep), '*'+seach_word_N+'*.dat')))):
            x=[]
            for f in filelist:
                if fnmatch.fnmatch(f, '*N'+str(i)+'.dat'):  
                    infile = data_path(dir+rep+f)
                    x = np.concatenate((x,np.genfromtxt(infile)))
                    
            
            (mean, s ,k) = block(x) #blocking analysis
            std = sqrt(s)
            data ={'N':[N_val[i]],'Mean':[mean], 'signa_block':[std[k]], 'std':[std[0]]}
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
    print("ATTENTION! The scripts analyzes all the scripts contained in the varying_alpha folder! If you just reduced the number of alphas you should manually delete files.")
    var_alpha(dir)
if selector==2: #Varyng dt
    print("ATTENTION! The scripts analyzes all the scripts contained in the varying_dt folder! If you just reduced the number of dt values you should manually delete files.")
    var_dt(dir)
if selector==3: #Varyng N
    print("ATTENTION! The scripts analyzes all the scripts contained in the varying_N folder! If you just reduced the number of N values you should manually delete files.")
    var_N(dir)


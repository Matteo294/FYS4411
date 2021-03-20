# Common imports
import os, fnmatch
import csv
import re
import multiprocessing as mp
import pandas as pd
from pandas import DataFrame
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
ap = argparse.ArgumentParser()
ap.add_argument("-p", "--Parallel_simulation", required=False, help=" 0--> standard_file_analysis; 1--> parallel_file_analysis")
ap.add_argument("-c", "--Program_selector", required=False, help="Chose Program:  \
    \n 0--> Simple case one file analysis;\
    \n 1--> multiple file varying alpha; \
    \n 2--> multiple file varying dt;\
    \n 3--> multiple file varying N")
ap.add_argument("-s", "--Save_figure", required=False, help="to save images in a default repository write 'save' or 'y' ")
args = vars(ap.parse_args())

if args['Parallel_simulation']: 
    dir = './Data/parallel'
else:
    dir = "./Data/standard"

if args['Program_selector'] == None: #default value (basic analysis from 1 file, and it doesn't save the figure)
    selector = 0
    fig = False #don't save
else :
    selector =int(args['Program_selector'])
    fig = False
    if args['Save_figure'] is not None or args['Save_figure']=="n":
        fig = True #save figure

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
    # preliminaries
    n = len(x)
    d = int(log2(n))
    s, gamma = zeros(d), zeros(d)
    #print(x)
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

    (mean, var, k) = block(x) 
    std = sqrt(var)
    print("std[0]= ", std[0] )
    data ={'Mean':[mean], 'STDev':[std[k]]}
    #frame = pd.DataFrame(data,index=['Values'])
    print(data)
    print ("\n=========================================\n")
    plt.plot(arange(0, len(std), 1), std)
    plt.grid()
    plt.show()
        #if fig == True:
         #   savefigure(dname,std,"0")
        #plt.show()
        #if fig == True:
         #   savefigure(dname,std,"0")
        #plt.show()
        #if fig == True:
         #   savefigure(dname,std,"0")


#Varying alpha 
def var_alpha(dir):
    rep="/varying_alpha/" 
    filelist =sorted_alphanumeric(fnmatch.filter(os.listdir(dir+rep), "*.dat"))

    with open(os.path.join(dir+rep, "post_analysis_alpha"+'.csv'), "w", newline='') as fcsv:
        writer = csv.writer(fcsv,delimiter =',')
        writer.writerow(["energy","std"])
        for i in range(0,int(len(fnmatch.filter(os.listdir(dir+rep), '*core0*.dat')))): # the variable i cycles on the different alpha
            x=[]
            for f in filelist:   
                if fnmatch.fnmatch(f, '*alpha'+str(i)+'.dat'): # here we collect all the data from file that represent the same alpha 
                    infile = data_path(dir+rep+f)
                    x = np.concatenate((x,np.genfromtxt(infile)))
            
            (mean, var ,k) = block(x) #blocking analysis
            std = sqrt(var)
            data ={'Mean':[mean], 'STDev':[std[k]]}
            writer = csv.writer(fcsv,delimiter =',')
            writer.writerow([mean,std[k]])
            #frame = pd.DataFrame(data,index=['Values'])
            print("alpha nr.", i, " --> ", data)

        if fig == True:
            savefigure(dir +"/Figures"+rep,std,"alpha of data set " + str(i))
        #Decomment these row below to have a fast view of the charts     
        #plt.plot(arange(0, len(std), 1), std)
        #plt.show(block=False)
        #plt.pause(2)
        #plt.close()
    fcsv.close()

def var_dt(dir):
    rep= "/varying_dt/"
    filelist=sorted_alphanumeric(os.listdir(dir+rep))
    
    with open(os.path.join(dir+rep, "post_analysis_dt"+'.csv'), "w", newline='') as fcsv:
        writer = csv.writer(fcsv,delimiter =',')
        writer.writerow(["energy","STD"])
        for i in range(0,int(len(fnmatch.filter(os.listdir(dir+rep), '*core0*.dat')))):
            x=[]
            for f in filelist:
                if fnmatch.fnmatch(f, '*dt'+str(i)+'.dat'):  
                    infile = data_path(dir+rep+f)
                    x = np.concatenate((x,np.genfromtxt(infile)))
            
            (mean, var ,k) = block(x) 
            std = sqrt(var)
            data ={'Mean':[mean], 'STDev':[std[k]]}
            writer = csv.writer(fcsv,delimiter =',')
            writer.writerow([mean,std[k]])
            #frame = pd.DataFrame(data,index=['Values'])
            print("dt nr.", i, " --> ", data)
            if fig == True:
                 savefigure(dir+"/Figures"+rep,std,f)
        #Decomment these row below to have a fast view of the charts 
        #plt.plot(arange(0, len(std), 1), std)
        #plt.show(block=False)
        #plt.pause(3)
        #plt.close()

def var_N(dir):
    rep= "/varying_N/"
    filelist =sorted_alphanumeric(os.listdir(dir+rep))

    with open(os.path.join(dir+rep, "post_analysis_N"+'.csv'), "w", newline='') as fcsv:
        writer = csv.writer(fcsv,delimiter =',')
        writer.writerow(["energy","STD"])
        for i in range(0,int(len(fnmatch.filter(os.listdir(dir+rep), '*core0*.dat')))):
            x=[]
            for f in filelist:
                if fnmatch.fnmatch(f, '*N'+str(i)+'.dat'):  
                    infile = data_path(dir+rep+f)
                    x = np.concatenate((x,np.genfromtxt(infile)))
            (mean, var ,k) = block(x) 
            std = sqrt(var)
            data ={'Mean':[mean], 'STDev':[std[k]]}
            writer = csv.writer(fcsv,delimiter =',')
            writer.writerow([mean,std[k]])
            #frame = pd.DataFrame(data,index=['Values'])
            print("N nr.", i, " --> ", data)
            if fig == True:
                   savefigure(dir+"/Figures"+rep ,std,f)
        #Decomment these row below to have a fast view of the charts 
        #plt.plot(arange(0, len(std), 1), std)
        #plt.show(block=False)
        #plt.pause(3)
        #plt.close()

#sort alphanurecaly the files
def sorted_alphanumeric(data):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(data, key=alphanum_key)
#save figure in the right repository
def savefigure(dirname,std,f): 
    font = {'fontname':'serif'}
    plt.plot(arange(0, len(std), 1), std)
    plt.grid()
    plt.ylabel('Standard error ', **font)
    plt.xlabel('Number of Blocking Transformations ', **font)
    plt.title('Standard error estimation by blocking method '+ f, **font)
    #plt.savefig(dirname + "/blocking"+f[-1]+".eps")
    print(f)
    plt.savefig(dirname + "/blocking"+str(*re.findall(r'\d+',f))+".png")
    plt.clf()

################################################################
#################################################################
#MAIN
#################################################################
#Case 
if selector==0: #Simple run
    simplerun(dir)
if selector==1: #Varyng alpha
    var_alpha(dir)
if selector==2: #Varyng dt
    var_dt(dir)
if selector==3: #Varyng N
    var_N(dir)

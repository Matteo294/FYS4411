# Common imports
import os
import re
import pandas as pd
from pandas import DataFrame
import matplotlib.pyplot as plt
import numpy as np
import sys
import argparse
import time #to implement the running time 
from numpy import log2, zeros, mean, var, sum, loadtxt, arange, array, cumsum, dot, transpose, diagonal, sqrt
from numpy.linalg import inv

#instruction how to run the program 
ap = argparse.ArgumentParser()
ap.add_argument("-c", "--Program_selector", required=False, help="Chose Program:  \
    \n 0--> Simple case one file analysis;\
    \n 1--> multiple file varying alpha; \
    \n 2--> multiple file varying N")
ap.add_argument("-s", "--Save_figure", required=False, help="to save images in a default repository write 'save' or 'y' ")
args = vars(ap.parse_args())


if args['Program_selector'] == None: #default value (basic analysis from 1 file, and it doesn't save the figure)
    selector = 0
    fig = False #don't save
else :
    selector =int(args['Program_selector'])
    fig = False
    if args['Save_figure'] is not None or "n":
        fig = True #save figure

#args = vars(ap.parse_args())
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# Choose the repository 
DATA_ID = "./"
def data_path(dat_id):
    return os.path.join(DATA_ID, dat_id)

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
    
    return mu, s[0:k]
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
    plt.savefig(dirname + "/blocking"+f[-1]+".png")
    plt.clf()

#################################################################
#Case 

if selector==0 : #Simple (default) case 1 file "energyateverystep"
    print ("\n=========================================\n")
    print("File:  energyateverystep.dat")
    infile = data_path(dname+"/energyateverystep.dat")
    #os.chdir("var_alpha")
    x = np.genfromtxt(infile)
    (mean, var) = block(x) 
    std = sqrt(var)
    data ={'Mean':[mean], 'STDev':[std[-1]]}
    frame = pd.DataFrame(data,index=['Values'])
    print(frame)
    if fig == True:
        savefigure(dname,std,"0")
#Case Varyng alpha
if selector==1 :  
    filelist =sorted_alphanumeric(os.listdir("./var_alpha"))
    st=0
    for f in filelist:
        print ("\n=========================================\n")
        print("File:", f)
        infile = data_path(dname+"/var_alpha/"+f)
        #os.chdir("var_alpha")
        x = np.genfromtxt(infile)
        (mean, var) = block(x) 
        std = sqrt(var)

        data ={'Mean':[mean], 'STDev':[std[-1]]}
        frame = pd.DataFrame(data,index=['Values'])
        print(frame)
        if fig == True:
            savefigure(dname +"/images_var_alpha/",std,"alpha of data set" + str(st))
            st+=1 
        #Decomment these row below to have a fast view of the charts     
        plt.plot(arange(0, len(std), 1), std)
        plt.show(block=False)
        plt.pause(3)
        plt.close()

#Varyng N
if selector==2: 
    filelist =sorted_alphanumeric(os.listdir("./var_N"))
    filelist.sort()
    for f in filelist:
        print ("\n=========================================\n")
        print("File:", f)
        infile = data_path(dname+"/var_N/"+f)
        x = np.genfromtxt(infile)    
        (mean, var) = block(x) 
        std = sqrt(var)
        data ={'Mean':[mean], 'STDev':[std[-1]]}
        frame = pd.DataFrame(data,index=['Values'])
        print(frame)
        if fig == True:
            savefigure(dname+"/images_var_N/",std,f)
        #Decomment these row below to have a fast view of the charts 
        #plt.plot(arange(0, len(std), 1), std)
        #plt.show(block=False)
        #plt.pause(3)
        #plt.close()
print ("\n=========================================\n")

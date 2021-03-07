# Common imports
import os
from numpy import std, mean, concatenate, arange, loadtxt, zeros, ceil 
from numpy.random import randint
from time import time
# Where to save the figures and data files
DATA_ID = "./" 
def data_path(dat_id):
        return os.path.join(DATA_ID, dat_id) 
#print(os.listdir())
infile = open(data_path("energy.dat"),'r')

def tsboot(data,statistic,R,l):
        t = zeros(R);
        n = len(data); 
        k = int(ceil(float(n)/l)); 
        inds = arange(n); 
        t0 = time();
        # time series bootstrap
        print(ceil(float(n)/l))
        for i in range(R):

                # construct bootstrap sample from
                # k chunks of data. The chunksize is l
                _data = data[randint(0,n,n)]; 
                t[i] = statistic(_data)
                #print(t[i])
                #print(_data)
        #print(_data)
# analysis
        print ("Runtime: %g sec" % (time()-t0)) 
        print ("Bootstrap Statistics :") 
        print ("original bias std. error")
        #print(std(t))
        #print(mean(t))
        print(statistic(data))
        print ("%8g %14g %15g" % (statistic(data), mean(t) - statistic(data), std(t)))
        return t
# Read in data
X = loadtxt(infile)
# statistic to be estimated. Takes two args. # arg1: the data
def stat(data):
        return mean(data)
t = tsboot(X, stat, 100, 1)
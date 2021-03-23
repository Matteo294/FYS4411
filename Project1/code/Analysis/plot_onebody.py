import pandas as pd
import matplotlib 
import matplotlib.pyplot as plt

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'


#data0 = pd.read_csv(r'./Data/parallel/onebody_density/counts0.csv', usecols = ["r", "counts"])
#data1 = pd.read_csv(r'./Data/parallel/onebody_density/counts1.csv', usecols = ["r", "counts"])
#data2 = pd.read_csv(r'./Data/parallel/onebody_density/counts2.csv', usecols = ["r", "counts"])
#data3 = pd.read_csv(r'./Data/parallel/onebody_density/counts3.csv', usecols = ["r", "counts"])

data00 = pd.read_csv(r'./Data/parallel/onebody_density/counts00.csv', usecols = ["r", "counts"])
data01 = pd.read_csv(r'./Data/parallel/onebody_density/counts01.csv', usecols = ["r", "counts"])
data02 = pd.read_csv(r'./Data/parallel/onebody_density/counts02.csv', usecols = ["r", "counts"])
data03 = pd.read_csv(r'./Data/parallel/onebody_density/counts03.csv', usecols = ["r", "counts"])

data10 = pd.read_csv(r'./Data/parallel/onebody_density/counts10.csv', usecols = ["r", "counts"])
data11 = pd.read_csv(r'./Data/parallel/onebody_density/counts11.csv', usecols = ["r", "counts"])
data12 = pd.read_csv(r'./Data/parallel/onebody_density/counts12.csv', usecols = ["r", "counts"])
data13 = pd.read_csv(r'./Data/parallel/onebody_density/counts13.csv', usecols = ["r", "counts"])

data20 = pd.read_csv(r'./Data/parallel/onebody_density/counts20.csv', usecols = ["r", "counts"])
data21 = pd.read_csv(r'./Data/parallel/onebody_density/counts21.csv', usecols = ["r", "counts"])
data22 = pd.read_csv(r'./Data/parallel/onebody_density/counts22.csv', usecols = ["r", "counts"])
data23 = pd.read_csv(r'./Data/parallel/onebody_density/counts23.csv', usecols = ["r", "counts"])

data30 = pd.read_csv(r'./Data/parallel/onebody_density/counts30.csv', usecols = ["r", "counts"])
data31 = pd.read_csv(r'./Data/parallel/onebody_density/counts31.csv', usecols = ["r", "counts"])
data32 = pd.read_csv(r'./Data/parallel/onebody_density/counts32.csv', usecols = ["r", "counts"])
data33 = pd.read_csv(r'./Data/parallel/onebody_density/counts33.csv', usecols = ["r", "counts"])

data0 = data00/4 + data01/4 + data02/4 + data03/4
data1 = data10/4 + data11/4 + data12/4 + data13/4
data2 = data20/4 + data21/4 + data22/4 + data23/4
data3 = data30/4 + data31/4 + data32/4 + data33/4


plt.figure(figsize=(10,8))
alpha=0.5
for i in range(len(data0.r)):
    if (i==0):
        ##plt.plot([data0.r[i], data0.r[i]], [0, data0.counts[i]], color='blue', label=r'$metropolis$')
        plt.plot([data0.r[i], data0.r[i]], [0, data0.counts[i]], alpha=alpha, color='red', label=r'$a=0.0000$')
        plt.plot([data1.r[i], data1.r[i]], [0, data1.counts[i]], alpha=0.4, color='yellow', label=r'$a=0.0043$')
        plt.plot([data2.r[i], data2.r[i]], [0, data2.counts[i]], alpha=0.7, color='deepskyblue', label=r'$a=0.043$')
        plt.plot([data3.r[i], data3.r[i]], [0, data3.counts[i]], alpha=alpha, color='green', label=r'$a=0.43$')
    else :
        #plt.plot([data0.r[i], data0.r[i]], [0, data0.counts[i]], color='blue')
        plt.plot([data0.r[i], data0.r[i]], [0, data0.counts[i]], alpha=alpha, color='red')
        plt.plot([data1.r[i], data1.r[i]], [0, data1.counts[i]], alpha=0.4, color='yellow')
        plt.plot([data2.r[i], data2.r[i]], [0, data2.counts[i]], alpha=0.7, color='deepskyblue')
        plt.plot([data3.r[i], data3.r[i]], [0, data3.counts[i]], alpha=alpha, color='green')



plt.xlabel(r'$r$', fontsize=22, labelpad=15)
plt.ylabel('Normalized counts', fontsize=22, labelpad=15)
plt.legend(fontsize=16)
ax = plt.gca()
ax.tick_params(axis='both', which='major', pad=5, labelsize=16)
plt.grid()
plt.show()

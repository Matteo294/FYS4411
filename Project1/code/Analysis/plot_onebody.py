import pandas as pd
import matplotlib 
import matplotlib.pyplot as plt

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'


data0 = pd.read_csv(r'./Data/parallel/onebody_density/counts0.csv', usecols = ["r", "counts"])
data1 = pd.read_csv(r'./Data/parallel/onebody_density/counts1.csv', usecols = ["r", "counts"])
data2 = pd.read_csv(r'./Data/parallel/onebody_density/counts2.csv', usecols = ["r", "counts"])
data3 = pd.read_csv(r'./Data/parallel/onebody_density/counts3.csv', usecols = ["r", "counts"])

data00 = pd.read_csv(r'./Data/parallel/onebody_density/counts00.csv', usecols = ["r", "counts"])
data01 = pd.read_csv(r'./Data/parallel/onebody_density/counts01.csv', usecols = ["r", "counts"])
data02 = pd.read_csv(r'./Data/parallel/onebody_density/counts02.csv', usecols = ["r", "counts"])
data03 = pd.read_csv(r'./Data/parallel/onebody_density/counts03.csv', usecols = ["r", "counts"])

data = data0/4 + data1/4 + data2/4 + data3/4
data0 = data00/4 + data01/4 + data02/4 + data03/4


plt.figure(figsize=(10,8))
for i in range(len(data.r)):
    if (i==0):
        plt.plot([data0.r[i], data0.r[i]], [0, data0.counts[i]], color='blue', label=r'$metropolis$')
        plt.plot([data.r[i], data.r[i]], [0, data.counts[i]], color='red', label=r'$importance$')
    else :
        plt.plot([data0.r[i], data0.r[i]], [0, data0.counts[i]], color='blue')
        plt.plot([data.r[i], data.r[i]], [0, data.counts[i]], color='red')



plt.xlabel(r'$r$', fontsize=22, labelpad=15)
plt.ylabel('Normalized counts', fontsize=22, labelpad=15)
plt.legend(fontsize=16)
ax = plt.gca()
ax.tick_params(axis='both', which='major', pad=5, labelsize=16)
plt.grid()
plt.show()

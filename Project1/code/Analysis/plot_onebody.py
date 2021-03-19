import pandas as pd
import matplotlib 
import matplotlib.pyplot as plt

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'


data = pd.read_csv(r'./data/onebody_density.csv', usecols = ["r", "counts"])



plt.figure(figsize=(10,8))
for i in range(len(data.r)):
    if (i==0):
        plt.plot([data.r[i], data.r[i]], [0, data.counts[i]], color='blue', label=r'$a=0.0043$')
    else :
        plt.plot([data.r[i], data.r[i]], [0, data.counts[i]], color='blue')



plt.xlabel(r'$r$', fontsize=22, labelpad=15)
plt.ylabel('Normalized counts', fontsize=22, labelpad=15)
plt.legend(fontsize=16)
ax = plt.gca()
ax.tick_params(axis='both', which='major', pad=5, labelsize=16)
plt.grid()
plt.show()

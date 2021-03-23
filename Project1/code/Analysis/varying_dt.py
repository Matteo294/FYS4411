import numpy as np 
from matplotlib import pyplot as plt 
import pandas as pd
import matplotlib

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

x = pd.read_csv('./Data/parallel/varying_dt/varying_dt.csv')
y = pd.read_csv('./Data/parallel/varying_dt/post_analysis_dt.csv')
data = pd.concat([x,y.reindex(x.index)], axis=1 ) # I concatenate the data in a unique dataframe


plt.figure(figsize=(10,8))
plt.errorbar(data['dt'], data['energy'], yerr=data['STD'], ls=None, marker='.', color = 'blue')
#plt.plot(xline, exactval, linewidth=0.7, color='blue', alpha=0.8)
#plt.scatter(data['dt'], data['energy'], color='red', s=12, zorder=2.5)
#plt.errorbar(data['dt'], data['energy'], yerr=data['std'], ls='none' , color='red', zorder=2.5)
plt.xlabel(r'$\delta t$', fontsize=22, labelpad=15)
plt.ylabel(r'Energy [$\hbar \omega_{ho}$]', fontsize=22, labelpad=15)
plt.legend(fontsize=16)
plt.yscale('log')
ax = plt.gca()
ax.tick_params(axis='both', which='major', pad=5, labelsize=16)
plt.grid()
#plt.savefig('./Figures/energy_varying_dt.eps')
plt.show()

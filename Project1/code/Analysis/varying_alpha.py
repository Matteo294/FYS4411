import numpy as np 
from matplotlib import pyplot as plt 
import pandas as pd
import matplotlib

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

x = pd.read_csv('./Data/parallel/varying_alpha/varying_alpha.csv')
y = pd.read_csv('./Data/parallel/varying_alpha/post_analysis_alpha.csv')
data = pd.concat([x,y.reindex(x.index)], axis=1 ) # I concatenate the data in a unique dataframe
xline = np.linspace(0.3, 0.7, 1000)
exactval = [30*(x/2 + 1/8/x) for x in xline]


plt.figure(figsize=(10,8))
plt.plot(xline, exactval, linewidth=1.8, color='mediumspringgreen', alpha=0.8, label='Model')
plt.scatter(data['alpha'], data['energy'], color='red', s=12, label='Importance', zorder=2.5)
plt.errorbar(data['alpha'], data['energy'], yerr=data['std'], ls='none' , color='red', zorder=2.5)
plt.xlabel(r'$\alpha$', fontsize=22, labelpad=15)
plt.ylabel('Energy', fontsize=22, labelpad=15)
ax = plt.gca()
plt.legend(fontsize=16)
ax.tick_params(axis='both', which='major', pad=5, labelsize=16)
plt.grid()
plt.savefig('./Figures/varying_alpha_noninteract_importance.eps')
plt.show()

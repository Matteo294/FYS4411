import numpy as np 
from matplotlib import pyplot as plt 
import pandas as pd
import matplotlib

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

alpha = pd.read_csv('./Data/parallel/varying_alpha/varying_alpha.csv')
energy, std = pd.read_csv('./Data/parallel/varying_alpha/post_analysis.csv')
#xline = np.linspace(0.3, 0.7, 1000)
#exactval = [30*(x/2 + 1/8/x) for x in xline]

plt.figure(figsize=(10,8))
#plt.plot(xline, exactval, linewidth=1.8, color='mediumspringgreen', alpha=0.8, label='Model')
plt.errorbar(alpha, energy, ls='none', yerr=std, color='red', zorder=2.5)
plt.scatter(alpha, energy, color='red', s=12, label='Metropolis', zorder=2.5)
plt.xlabel(r'$\alpha$', fontsize=22, labelpad=15)
plt.ylabel('Energy', fontsize=22, labelpad=15)
ax = plt.gca()
plt.legend(fontsize=16)
ax.tick_params(axis='both', which='major', pad=5, labelsize=16)
plt.grid()
#plt.savefig('./figures/varying_alpha_noninteract.eps')
plt.show()
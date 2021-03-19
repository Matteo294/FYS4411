import numpy as np 
from matplotlib import pyplot as plt 
import pandas as pd
import matplotlib

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

data = pd.read_csv('./data/varying_alpha.csv')

# Non interacting, Metropolis
alphas = np.arange(0.3, 0.71, 0.02)
val = [17.0216, 16.5288, 16.1293, 15.8140, 15.5648, 15.3676, 15.2408, 15.1185, 15.0511, 15.0119, 15.0, 15.0118, 15.0433, 15.0949, 15.1686, 15.2436, 15.3537, 15.4424, 15.6039, 15.7269, 15.8442]
std = [0.0165,   0.0139,  0.0118,  0.0098,  0.0080,  0.0064,  0.0053,  0.0037,  0.0023,  0.0011,  0.0,  0.0011,  0.0021,  0.0033,  0.0041,  0.0052,  0.0061,  0.0070,  0.0071,  0.0078,  0.0095]
xline = np.linspace(0.3, 0.7, 1000)
exactval = [30*(x/2 + 1/8/x) for x in xline]

plt.figure(figsize=(10,8))
plt.plot(xline, exactval, linewidth=1.8, color='mediumspringgreen', alpha=0.8, label='Model')
plt.errorbar(alphas, val, ls='none', yerr=std, color='red', zorder=2.5)
plt.scatter(alphas, val, color='red', s=12, label='Metropolis', zorder=2.5)
plt.xlabel(r'$\alpha$', fontsize=22, labelpad=15)
plt.ylabel('Energy', fontsize=22, labelpad=15)
ax = plt.gca()
plt.legend(fontsize=16)
ax.tick_params(axis='both', which='major', pad=5, labelsize=16)
plt.grid()
#plt.savefig('./figures/varying_alpha_noninteract.eps')
plt.show()
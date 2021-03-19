import numpy as np 
from matplotlib import pyplot as plt 
import pandas as pd
import matplotlib

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

data = pd.read_csv('./data/varying_dt.csv')

plt.figure(figsize=(10,8))
plt.plot(data['dt'], data['acceptance'], linewidth=2.2, color='cornflowerblue')
plt.xscale('log')
xlab = plt.xlabel(r'$\delta t$', fontsize=22, labelpad=15)
ylab = plt.ylabel('Acceptance', fontsize=22, labelpad=15)
ax = plt.gca()
ax.tick_params(axis='both', which='major', pad=5, labelsize=16)
plt.grid()
plt.savefig('./figures/varying_dt.eps')
plt.show()

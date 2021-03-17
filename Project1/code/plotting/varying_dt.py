import numpy as np 
from matplotlib import pyplot as plt 
import pandas as pd
import matplotlib

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

data = pd.read_csv('./data/varying_dt.csv')

plt.plot(data['dt'], data['acceptance'], linewidth=1.8, color='cornflowerblue')
plt.xscale('log')
xlab = plt.xlabel(r'$\delta t$', fontsize=16, labelpad=15)
ylab = plt.ylabel('Acceptance', fontsize=16, labelpad=15)
ax = plt.gca()
ax.tick_params(axis='both', which='major', pad=5)
#plt.grid()
plt.show()
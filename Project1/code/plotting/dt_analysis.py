import numpy as np 
from matplotlib import pyplot as plt 
import pandas as pd
import matplotlib

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

data = pd.read_csv('./data/varying_dt.csv')

plt.plot(data['dt'], data['energy'], linewidth=1.8, color='black')
plt.xlabel(r'dt', fontsize=14)
plt.ylabel('E / N', fontsize=14)
plt.grid()
plt.show()
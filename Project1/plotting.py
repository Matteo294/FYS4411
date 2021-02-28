import pandas as pd 
from matplotlib import pyplot as plt 

data = pd.read_csv('results.csv')

plt.plot(data['alpha'], data['energy'], linewidth=1.8)

data2 = pd.read_csv('resultfile.csv')

plt.plot(data2['alpha'], data2['averEnergy'])

plt.show()
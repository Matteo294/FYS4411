import pandas as pd 
from matplotlib import pyplot as plt 

data = pd.read_csv('results.csv')

plt.plot(data['alpha'], data['energy'], linewidth=1.8)
plt.show()
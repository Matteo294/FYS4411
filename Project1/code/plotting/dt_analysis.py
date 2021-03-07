import pandas as pd 
from matplotlib import pyplot as plt 

data = pd.read_csv('./data/dt_analysis.csv')

plt.plot(data['dt'], data['EL'], linewidth=1.8, color='black')
plt.xlabel('dt', fontsize=14)
plt.ylabel('Local Energy', fontsize=14)
plt.show()
import pandas as pd
import numpy as np
import matplotlib 
import matplotlib.pyplot as plt

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'


data = np.genfromtxt(r'./Data/standard/singlerun/energyateverystep.dat')

print(data)
plt.figure(figsize=(10,8))
plt.plot(range(len(data)), data)
print(min(data))

'''
plt.xlabel(r'$r$', fontsize=22, labelpad=15)
plt.ylabel('Normalized counts', fontsize=22, labelpad=15)
plt.legend(fontsize=16)
ax = plt.gca()
ax.tick_params(axis='both', which='major', pad=5, labelsize=16)
'''
plt.grid()
plt.show()

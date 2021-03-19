import pandas as pd
import matplotlib 
import matplotlib.pyplot as plt

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

plt.figure(figsize=(10,8))
colors = ['blue', 'red', 'yellow', 'green']
alphas = [0.5, 0.5, 0.5, 0.5]
labels = [r'$a=0.0$', r'$a=0.0043$', r'$a=0.043$', r'$a=0.43$']

for j in range(4):
    data = pd.read_csv(r'./data/onebody_density' + str(j) + '.csv', usecols = ["r", "counts"])

    for i in range(len(data.r)):
        if (i==0):
            plt.plot([data.r[i], data.r[i]], [0, data.counts[i]], color=colors[j], label=labels[j], alpha=alphas[j])
        else :
            plt.plot([data.r[i], data.r[i]], [0, data.counts[i]], color=colors[j], alpha=alphas[j])



plt.xlabel(r'$r$', fontsize=22, labelpad=15)
plt.ylabel('Normalized counts', fontsize=22, labelpad=15)
plt.legend(fontsize=16)
ax = plt.gca()
ax.tick_params(axis='both', which='major', pad=5, labelsize=16)
plt.grid()
plt.savefig('./figures/onebody_density.eps')
plt.show()


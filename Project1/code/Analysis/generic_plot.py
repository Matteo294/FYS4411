import pandas as pd
import numpy as np
import matplotlib 
import matplotlib.pyplot as plt

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

'''
y_axis = []
y_axis2 = []
for i in range(4):
    y_axis = np.concatenate((y_axis, np.genfromtxt('./Data/parallel/singlerun/energyateverystep' + str(i) +'_tenere.dat')))
    y_axis2 = np.concatenate((y_axis2, np.genfromtxt('./Data/parallel/singlerun/energyateverystep1' + str(i) +'_tenere.dat')))

y_axis=y_axis[0:2**13]
y_axis2=y_axis2[0:2**13]
x_axis = range(len(y_axis))

plt.figure(figsize=(10,8))
plt.plot(x_axis, y_axis2, linewidth = 1.0, color='red', label=r'$\delta t = 0.1$')
plt.plot(x_axis, y_axis, linewidth = 1.0, color='blue', label=r'$\delta t = 0.001$')
plt.xlabel(r'$N_{steps}$', fontsize=22, labelpad=15)
plt.ylabel(r'$E_L$ [$\hbar \omega_{ho}$]', fontsize=22, labelpad=15)
plt.legend(fontsize=16)
ax = plt.gca()
ax.tick_params(axis='both', which='major', pad=5, labelsize=16)
plt.grid()
#plt.savefig('./Figures/correlation_different_deltat.eps')
plt.show()
'''





#PLOT FOR ENERGY DIVIDED BY NUMBER OF PARTICLES IN THE INTERACTING CASE N=10, 50, 100

x_axis = [0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7]



energy10 = [27.588704850893023, 25.90947811525536, 24.969203620384707, 24.52252123904419,24.399054127834315, 24.525824823119176,
24.821768408590316, 25.265636166320817, 25.834827534567363]
energy10 = [x/10 for x in energy10]
denergy10 = [0.022856132618976167, 0.015595915826757149, 0.008418565621906742, 0.0040289201802913704, 0.00030743248558163453, 
0.003557179588239271, 0.006425119214188396, 0.008672248758691672, 0.01104995550586971]
denergy10 = [x/10 for x in denergy10]

'''
energy50 = [142.8, 134.53, 129.81, 127.71, 127.287, 129.90, 135.54]
energy50 = [x/50 for x in energy50]
denergy50 = [0.1, 0.09, 0.04, 0.02, 0.006, 0.03, 0.06]
denergy50 = [x/50 for x in denergy50]

energy100 = [296.8, 271.2, 266.51, 272.8, 285.1]
energy100 = [x/100 for x in energy100]
denergy100 = [0.4, 0.1, 0.07, 0.2, 0.2]
denergy100 = [x/100 for x in denergy100]
'''

plt.figure(figsize=(10,8))
plt.errorbar(x_axis, energy10, yerr=denergy10, marker='.', markersize=12, ls='-' , color='red', label='N=10')
#plt.errorbar(x_axis, energy50, yerr=denergy50, marker='.', markersize=12, ls='-' , color='blue', label='N=50')
#plt.errorbar(x_axis, energy100, yerr=denergy100, marker='.', markersize=12, ls='-' , color='purple', label='N=100')

plt.xlabel(r'$\alpha$', fontsize=22, labelpad=15)
plt.ylabel(r'Energy [$\hbar \omega_{ho}$]', fontsize=22, labelpad=15)
plt.legend(fontsize=16)
ax = plt.gca()
ax.tick_params(axis='both', which='major', pad=5, labelsize=16)
plt.grid()
plt.savefig('./Figures/energy_over_N.eps')
plt.show()




'''
x_axis = [0.3, 0.4, 0.5, 0.6, 0.7]

energy10 = [27.62, 24.986, 24.3992, 24.827, 25.84]
energy10 = [x/10 for x in energy10]
denergy10 = [0.02, 0.008, 0.0003, 0.007, 0.01]
denergy10 = [x/10 for x in denergy10]

energy50 = [142.8, 129.81, 127.287, 129.90, 135.54]
energy50 = [x/50 for x in energy50]
denergy50 = [0.1, 0.04, 0.006, 0.03, 0.06]
denergy50 = [x/50 for x in denergy50]

energy100 = [296.8, 271.2, 266.51, 272.8, 285.1]
energy100 = [x/100 for x in energy100]
denergy100 = [0.4, 0.1, 0.07, 0.2, 0.2]
denergy100 = [x/100 for x in denergy100]


plt.figure(figsize=(10,8))
plt.errorbar(x_axis, energy10, yerr=denergy10, marker='.', markersize=12, ls='-' , color='red', label='N=10')
plt.errorbar(x_axis, energy50, yerr=denergy50, marker='.', markersize=12, ls='-' , color='blue', label='N=50')
plt.errorbar(x_axis, energy100, yerr=denergy100, marker='.', markersize=12, ls='-' , color='purple', label='N=100')

plt.xlabel(r'$\alpha$', fontsize=22, labelpad=15)
plt.ylabel(r'Energy [$\hbar \omega_{ho}$]', fontsize=22, labelpad=15)
plt.legend(fontsize=16)
ax = plt.gca()
ax.tick_params(axis='both', which='major', pad=5, labelsize=16)
plt.grid()
plt.savefig('./Figures/energy_over_N.eps')
plt.show()
'''



# PLOT PER GRAFICO DI COMPARAZIONE PER N=10 3D INTERAGENTE E NON
'''
y_axis = pd.read_csv(r'./Data/parallel/varying_alpha/post_analysis_alpha_N10int_tenerlo.csv')
x_axis = pd.read_csv(r'./Data/parallel/varying_alpha/varying_alpha.csv')

y_axis2 = pd.read_csv(r'./Data/parallel/varying_alpha/post_analysis_alpha_N10nonint_tenerlo.csv')

plt.figure(figsize=(10,8))
plt.plot(x_axis['alpha'], y_axis['energy'], marker='.', markersize=12, color='red', label='Interacting', zorder=2.5)
plt.errorbar(x_axis['alpha'], y_axis['energy'], yerr=y_axis['std'], ls='none' , color='red', zorder=2.5)
plt.plot(x_axis['alpha'], y_axis2['energy'],marker='.', markersize=12, color='blue', label='Non-interacting', zorder=2.5)
plt.errorbar(x_axis['alpha'], y_axis2['energy'], yerr=y_axis2['std'], ls='none' , color='blue', zorder=2.5)

plt.xlabel(r'$\alpha$', fontsize=22, labelpad=15)
plt.ylabel(r'Energy [$\hbar \omega_{ho}$]', fontsize=22, labelpad=15)
plt.legend(fontsize=16)
ax = plt.gca()
ax.tick_params(axis='both', which='major', pad=5, labelsize=16)
plt.grid()
#plt.savefig('./Figures/varying_alpha_comparison_interacting_and_not.eps')
plt.show()
'''

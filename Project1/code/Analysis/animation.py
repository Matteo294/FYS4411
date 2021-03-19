from matplotlib import pyplot as plt
import numpy as np
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import animation
import pandas as pd
from matplotlib import cm

fig = plt.figure()
ax = p3.Axes3D(fig)

R = 10
K = 1
scale = 1

_data = pd.read_csv('./data/particle.csv')
_data = _data.transpose()
data = _data.to_numpy()
N = data.shape[1]
for i in range(data.shape[0]):
    for j in range(data.shape[1]):
        data[i,j]
RMAX = data.max()

r = np.linspace(-K*RMAX, K*RMAX, 100)
x, y = np.meshgrid(r, r)
R = x**2 + y**2
z = np.exp(-0.5*R)*scale

def update(num, data, line):
    line.set_data(data[:2, :num])
    line.set_3d_properties(data[2, :num])

surf = ax.plot_surface(x, y, z, cmap='gray', alpha=0.5, linewidth=0, antialiased=False)
line, = ax.plot(data[0, 0:1], data[1, 0:1], data[2, 0:1])

# Setting the axes properties
ax.set_xlim3d([-K*RMAX, K*RMAX])
ax.set_xlabel('X')

ax.set_ylim3d([-K*RMAX, K*RMAX])
ax.set_ylabel('Y')

ax.set_zlim3d([0, RMAX])
ax.set_zlabel('Z')

ani = animation.FuncAnimation(fig, update, N, fargs=(data, line), interval=10000/N, blit=False)
#ani.save('matplot003.gif', writer='imagemagick')

ax._axis3don = False
plt.show()

import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv(r'./data/onebody_density.csv', usecols = ["r", "counts"])


for i in range(len(data.r)):
    plt.plot([data.r[i], data.r[i]], [0, data.counts[i]])

plt.show()
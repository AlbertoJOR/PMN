import pandas as pd
import matplotlib.pyplot as plt


plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True


headers = ['step', 'ISING1','ISING2']
cells= ['ISING1','ISING2']

df = pd.read_csv('ising.csv', usecols=headers,)

df.set_index('step')
df.loc[0:175,cells].plot()
plt.xlabel("step")
plt.ylabel("Energía Ising")
plt.grid()


plt.show()
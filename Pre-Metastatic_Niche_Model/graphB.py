import pandas as pd
import matplotlib.pyplot as plt


plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True


headers = ['step', 'BMDC','CTC','T_CD8','M1','N1','M2','N2','PD_L1']
cells= [ 'CTC','BMDC','T_CD8','M1','N1','M2','N2','PD_L1']

df = pd.read_csv('result_cell_abundanceB.csv', usecols=headers,)

df.set_index('step')
df.loc[0:175,cells].plot()
plt.xlabel("step")
plt.ylabel("Número de agentes")


plt.show()
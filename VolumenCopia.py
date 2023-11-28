import numpy as np
import matplotlib.pyplot as plt

tm = np.genfromtxt("VolumenCopia.dat", delimiter=";", skip_header=1)						#Importo la matriz de temperaturas almacenado como vector
tm = np.transpose(tm)

fig, ax = plt.subplots()								#Inicialización de una figura con matplotlib
#ax.plot(tm[0],tm[1]*1000, label='Exterior')
ax.plot(tm[0],tm[4]*1000000, label='Cáscara')
ax.plot(tm[0],tm[5]*1000000, label='Albúmina')
ax.plot(tm[0],tm[6]*1000000, label='Yema')
ax.plot(tm[0],tm[7]*1000000, label='Embrión')

ax.legend(loc="lower right", title="Región")
ax.set(xlabel='Tiempo (Horas)', ylabel='Volumen (cm^3)',
	title='Volumen de cada región vs tiempo')
ax.grid(True)


fig.savefig('VolumenCopia.jpg')


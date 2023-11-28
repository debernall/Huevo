import numpy as np
import matplotlib.pyplot as plt

tm = np.genfromtxt("TemperaturaMedia.dat", delimiter=";", skip_header=1)						#Importo la matriz de temperaturas almacenado como vector
tm = np.transpose(tm)

fig, ax = plt.subplots()								#Inicialización de una figura con matplotlib
ax.plot(tm[0],tm[1], label='Exterior')
ax.plot(tm[0],tm[2], label='Nido')
ax.plot(tm[0],tm[3], label='Parche incubación')
ax.plot(tm[0],tm[4], label='Cáscara')
ax.plot(tm[0],tm[5], label='Albúmina')
ax.plot(tm[0],tm[6], label='Yema')
ax.plot(tm[0],tm[7], label='Embrión')
#plt.ylim(0,55)
plt.xlim(0,tm[0][-1])
ax.legend(loc="upper right", title="Región")
ax.set(xlabel='Tiempo (Horas)', ylabel='Temperatura (°C)',
	title='Temperatura media de cada región vs tiempo')
ax.grid(True)


fig.savefig('TemperaturaMedia.jpg')


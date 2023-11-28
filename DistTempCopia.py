import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


t0 = np.genfromtxt("Matriz0Copia.dat")						#Importo la matriz de temperaturas almacenado como vector


n = np.cbrt(np.shape(t0)[0]).astype(int)					#Cálculo de la magnitud de cada vector de la matriz temperatura
t0 = t0.reshape((n,n,n))							#Reorganizo el vector en una matriz de nxnxn
ang1=10
ang2=20
size=0.9/n									#Tamaño de la representación de los nodos

m = 25.0									#Distancia máxima en cada eje desde el origen	
d = (2*m/(n)) 									#Delta, separación entre cada nodo
m1 = np.arange(-m,m,d)                         					#Arreglo de tamaño n que sirve de base para la construcción de las matrices x,y,z
x, y, z = np.meshgrid(m1,m1,m1, indexing='ij')  				#Construcción de una red de puntos que me va a permitir graficar en un espacio una variable en un espacio 3 Dimensional

fig = plt.figure()								#Inicialización de una figura con matplotlib

ax1 = fig.add_subplot(1, 1, 1, projection='3d')
scatter1 = ax1.scatter3D(x.flatten(), y.flatten(), z.flatten(), c=t0.flatten()*1, cmap='inferno', s=size, vmin=0, vmax=100)		
ax1.view_init(ang1,ang2)								
ax1.set_title("0% Tiempo")
				

pcbar1 = plt.colorbar(scatter1, orientation='vertical', pad=0.1, location='right', label='Temp °C', ax=ax1)

plt.subplots_adjust(hspace=0.5)

ax1.tick_params(axis='both', labelsize=7)


fig.suptitle("Temperatura para diferentes tiempos")
fig.savefig('DistTempCopia.jpg')


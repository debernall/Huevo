import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


t0 = np.genfromtxt("Matriz0.dat")						#Importo la matriz de temperaturas almacenado como vector
t1 = np.genfromtxt("Matriz33.dat")
t2 = np.genfromtxt("Matriz66.dat")
t3 = np.genfromtxt("Matriz100.dat")

n = np.cbrt(np.shape(t0)[0]).astype(int)					#Cálculo de la magnitud de cada vector de la matriz temperatura
t0 = t0.reshape((n,n,n))							#Reorganizo el vector en una matriz de nxnxn
t1 = t1.reshape((n,n,n))
t2 = t2.reshape((n,n,n))
t3 = t3.reshape((n,n,n))
ang1=10
ang2=20
size=0.005/n									#Tamaño de la representación de los nodos

m = 25.0									#Distancia máxima en cada eje desde el origen	
d = (2*m/(n)) 									#Delta, separación entre cada nodo
m1 = np.arange(-m,m,d)                         					#Arreglo de tamaño n que sirve de base para la construcción de las matrices x,y,z
x, y, z = np.meshgrid(m1,m1,m1, indexing='ij')  				#Construcción de una red de puntos que me va a permitir graficar en un espacio una variable en un espacio 3 Dimensional

fig = plt.figure()								#Inicialización de una figura con matplotlib

ax1 = fig.add_subplot(2, 2, 1, projection='3d')
scatter1 = ax1.scatter3D(x.flatten(), y.flatten(), z.flatten(), c=t0.flatten()*1, cmap='inferno', s=size, vmin=0, vmax=55)		
ax1.view_init(ang1,ang2)								
ax1.set_title("t=0% Tiempo")
				
ax2 = fig.add_subplot(2,2,2, projection='3d')					
scatter2 = ax2.scatter3D(x.flatten(), y.flatten(), z.flatten(), c=t1.flatten()*1, cmap='inferno', s=size, vmin=0, vmax=55)		
ax2.view_init(ang1,ang2)								
ax2.set_title("t=33% Tiempo")

ax3 = fig.add_subplot(2,2,3, projection='3d')					
scatter3 = ax3.scatter3D(x.flatten(), y.flatten(), z.flatten(), c=t2.flatten()*1, cmap='inferno', s=size, vmin=0, vmax=55)		
ax3.view_init(ang1,ang2)								
ax3.set_title("t=66% Tiempo")

ax4 = fig.add_subplot(2,2,4, projection='3d')					
scatter4 = ax4.scatter3D(x.flatten(), y.flatten(), z.flatten(), c=t3.flatten()*1, cmap='inferno', s=size, vmin=0, vmax=55)	
ax4.view_init(ang1,ang2)								
ax4.set_title("t=100% Tiempo")

pcbar1 = plt.colorbar(scatter1, orientation='vertical', pad=0.1, location='right', label='Temp °C', ax=ax1)
pcbar2 = plt.colorbar(scatter2, orientation='vertical', pad=0.1, location='right', label='Temp °C', ax=ax2)
pcbar3 = plt.colorbar(scatter3, orientation='vertical', pad=0.1, location='right', label='Temp °C', ax=ax3)
pcbar4 = plt.colorbar(scatter4, orientation='vertical', pad=0.1, location='right', label='Temp °C', ax=ax4)

plt.subplots_adjust(hspace=0.5)

ax1.tick_params(axis='both', labelsize=7)
ax2.tick_params(axis='both', labelsize=7)
ax3.tick_params(axis='both', labelsize=7)
ax4.tick_params(axis='both', labelsize=7)

fig.suptitle("Temperatura para diferentes tiempos")
plt.savefig('Temperatura.jpg')


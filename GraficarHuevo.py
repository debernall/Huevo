import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


t0 = np.genfromtxt("Matriz0.dat")						#Importo la matriz de temperaturas almacenado como vector
t10 = np.genfromtxt("Matriz10.dat")
t30 = np.genfromtxt("Matriz30.dat")
t55 = np.genfromtxt("Matriz55.dat")

n = np.cbrt(np.shape(t0)[0]).astype(int)					#Cálculo de la magnitud de cada vector de la matriz temperatura
t0 = t0.reshape((n,n,n))							#Reorganizo el vector en una matriz de nxnxn
t10 = t10.reshape((n,n,n))
t30 = t30.reshape((n,n,n))
t55 = t55.reshape((n,n,n))
ang1=45
ang2=45
size=0.01/n

m = 25.0									#Distancia máxima en cada eje desde el origen	
d = (2*m/(n)) 									#Delta, separación entre cada nodo
m1 = np.arange(-m,m,d)                         					#Arreglo de tamaño n que sirve de base para la construcción de las matrices x,y,z
x, y, z = np.meshgrid(m1,m1,m1, indexing='ij')  		#Construcción de una red de puntos que me va a permitir graficar en un espacio una variable en un espacio 3 Dimensional

fig = plt.figure("")								#Inicialización de una figura con matplotlib

ax1 = fig.add_subplot(2,2,1, projection='3d')					#Primer subplot
scatter1 = ax1.scatter3D(x.flatten(), y.flatten(), z.flatten(), c=t0.flatten()*-1, cmap='hot', s=size)		#Gráfica de Temperatura, s es el tamaño de los nodos
ax1.view_init(ang1,ang2)								#Ángulos iniciales de observación de la gráfica generada
ax1.set_title("t=0")
#fig.colorbar(scatter1, label='Temperatura')					#Barra lateral que indica la escala de temperatura en función del color presentado en la gráfica

ax2 = fig.add_subplot(2,2,2, projection='3d')					#Primer subplot
scatter2 = ax2.scatter3D(x.flatten(), y.flatten(), z.flatten(), c=t10.flatten()*-1, cmap='hot', s=size)		#Gráfica de Temperatura, s es el tamaño de los nodos
ax2.view_init(ang1,ang2)								#Ángulos iniciales de observación de la gráfica generada
ax2.set_title("t=10")
#fig.colorbar(scatter2, label='Temperatura')					#Barra lateral que indica la escala de temperatura en función del color presentado en la gráfica

ax3 = fig.add_subplot(2,2,3, projection='3d')					#Primer subplot
scatter3 = ax3.scatter3D(x.flatten(), y.flatten(), z.flatten(), c=t30.flatten()*-1, cmap='hot', s=size)		#Gráfica de Temperatura, s es el tamaño de los nodos
ax3.view_init(ang1,ang2)								#Ángulos iniciales de observación de la gráfica generada
ax3.set_title("t=30")
#fig.colorbar(scatter3, label='Temperatura')					#Barra lateral que indica la escala de temperatura en función del color presentado en la gráfica

ax4 = fig.add_subplot(2,2,4, projection='3d')					#Primer subplot
scatter4 = ax4.scatter3D(x.flatten(), y.flatten(), z.flatten(), c=t55.flatten()*-1, cmap='hot', s=size)		#Gráfica de Temperatura, s es el tamaño de los nodos
ax4.view_init(ang1,ang2)								#Ángulos iniciales de observación de la gráfica generada
ax4.set_title("t=55")
#fig.colorbar(scatter4, label='Temperatura')					#Barra lateral que indica la escala de temperatura en función del color presentado en la gráfica

print(np.sum(t0))
print(np.sum(t10))
print(np.sum(t30))
print(np.sum(t55))

plt.savefig('Temperatura.jpg')



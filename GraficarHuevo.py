import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


Temp = np.genfromtxt("Matriz.dat")						#Importo la matriz de temperaturas almacenado como vector
n = np.cbrt(np.shape(Temp)[0]).astype(int)					#Cálculo de la magnitud de cada vector de la matriz temperatura
Temp = Temp.reshape((n,n,n))							#Reorganizo el vector en una matriz de nxnxn

m = 25.0									#Distancia máxima en cada eje desde el origen	
d = (2*m/(n)) 									#Delta, separación entre cada nodo

m1 = np.arange(-m,m,d)                         					#Arreglo de tamaño n que sirve de base para la construcción de las matrices x,y,z
x, y, z = np.meshgrid(m1,m1,m1, indexing='ij')  		#Construcción de una red de puntos que me va a permitir graficar en un espacio una variable en un espacio 3 Dimensional

fig = plt.figure()								#Inicialización de una figura con matplotlib
ax = fig.add_subplot(111, projection='3d')					#Primer subplot
scatter = ax.scatter3D(x.flatten(), y.flatten(), z.flatten(), c=datos.flatten(), cmap='hot', s=0.1/n)		#Gráfica de Temperatura, s es el tamaño de los nodos
ax.view_init(45, 45)								#Ángulos iniciales de observación de la gráfica generada
fig.colorbar(scatter, label='Temperatura')					#Barra lateral que indica la escala de temperatura en función del color presentado en la gráfica
plt.show()									#CAMBIAR A IMPORTAR A JPG 




import numpy as np
import matplotlib.pyplot as plt

datos = np.genfromtxt("Matriz.dat")
datos = datos.reshape((50,50,50))

n = 50                                       			#500 es un buen número, consumo 250MB
m = 25.0								# m es el valor máximo y mínimo alcanzado por la matriz del espacio. En este caso 	
								# se define un espacio formado entre -25 y 25 en x,y,z. Unidades m x 10^-4 
d = (2*m/(n)) 					# Delta, distancia entre cada punto. d=5. Uso np.int32 para reducir procesamiento

m1 = np.arange(-m,m,d)                         # m1 es un arreglo de tamaño n que sirve de base para la construcción de las matrices x,y,z
x, y, z = np.meshgrid(m1,m1,m1, indexing='ij')  

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
scatter = ax.scatter3D(x[:25,:,:], y[:25,:,:], z[:25,:,:], c=datos[:25,:,:], cmap='viridis', s=0.1)                 #Gráfica de la cáscara
ax.view_init(45, 45)
plt.show()




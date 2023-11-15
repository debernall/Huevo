import numpy as np
import matplotlib.pyplot as plt

datos = np.genfromtxt("Matriz.dat")
datos = datos.reshape((100,100,100))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
scatter = ax.voxels(datos, edgecolor='k')                 #Gráfica de la cáscara
ax.view_init(0, 90)
plt.show()



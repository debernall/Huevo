import numpy as np
import matplotlib.pylab as plt

# Malla cúbica que contiene las condiciones iniciales del sistema, en el interior se encontrará el huevo esférico.

# Inicialmente generaré una matriz de 3x3 con los valores 1 y 0 que me van a definir el espacio interno y externo al huevo.

# Tomo una sección transversal al huevo que me representa un corte transversal al plano x=0.
# La curva que me delimita esta sección cumple con la relación y^2+z^2=r^2
# Selecciono solamente la función y = sqrt(r**2-z**2)
# Dado que el huevo tendrá simetría de revolución, "y" me representará el radio para cada punto a lo largo de z
# Como me interesa usar un espacio mas grande 3D mas grande que el huevo, para los valores de z que no pertenecen a  

def radio(z,r):
    if np.abs(z)<r:
        radio = np.sqrt((r**2)-(z**2))
    else:
        radio = 0
    return radio

def espaciointerno(x,y,r):
    p = np.sqrt(x**2+y**2)
    if p<=r:
        out=0
    else:
        out=100
    return out
    
    
n = 50                                          #Tamaño de los vectores n=250
m = 25                                          #Determina el max y el min valor del espacio                                       
r = 24                                          #Radio max del huevo

m1 = np.linspace(m,-m,n)
m2 = np.linspace(-m,m,n)

y, x, z = np.meshgrid(m1,m1,m1)
      
M=np.ones((n,n,n))
for k in range(n):
    radio_z = radio(z[0][0][k],r)
    for j in range(n):
        for i in range(n):
            M[i][j][k] = espaciointerno(x[i][0][0],y[0][j][0],radio_z)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

scatter = ax.scatter3D(x, y, z, c=M, cmap='gray', s=0.001)

plt.show()
# for i in range(n):
#   for j in range(n):
#     if (X[i][j]**2+Y[i][j]**2)<r**2:
#       M[i][j] = 0


# Implementación del método de diferencias finitas para la ecuación de difusión de calor, usando una función

# Implementar la variación de la generación de calor en la yema

# Implementar el cambio de volumen de los componentes del huevo.


import numpy as np
import matplotlib.pylab as plt
import time

# Malla cúbica que contiene las condiciones iniciales del sistema, en el interior se encontrará el huevo esférico.

# Inicialmente generaré una matriz de 3x3 con los valores True y False que me van a definir el espacio interno y externo al huevo.

# Constantes de la simulación. 
# Construcción de la matriz cúbica del espacio, n es el tamaño de cada vector que compone la matriz del espacio.
# En este caso elijo n=50, por tanto n^3=125.000 son los puntos que componen la matriz del espacio
n = 500                                       #500 es un buen numero consumo 250MB

# m es el valor máximo y mínimo alcanzado por la matriz del espacio. En este caso se define un espacio formado entre
# -25 y 25 en x,y,z. Las unidades están dadas en micrómetros
m = 250                                          #Determina el max y el min valor del espacio   25000micrometros lo reduzco a 250
d = np.int(2*m/(n-1)) 

# En este caso cada punto representa un espacio equivalente a (xxxxxxx PENDIENTE POR CALCULAR)                                   

# Constantes para la formación del huevo, en este caso el huevo es esférico. r=24mm
R = 240                                          

# Creación de 3 matrices que representan las coordenadas x,y,z para cada punto dentro de la matriz de espacio
t0 = time.time()

m1 = np.arange(-m,m+1,d,dtype=np.int32)                        # m1 es un arreglo de tamaño n que sirve de base para la construcción de las matrices x,y,z
x, y, z = np.meshgrid(m1,m1,m1, indexing='ij')                 # Creación de las matrices x,y,z. El orden de los vectores se toma de manera conveniente

# Matriz del espacio      
# Creo una matriz similar a las máscaras de numpy usando lo siguiente:

# Tomo una sección transversal al huevo que me representa un corte transversal al plano x=0.
# La curva que me delimita esta sección cumple con la relación r^2+z^2=R^, donde r puede ser y o x o cualquier vector en el plano z.
# Selecciono solamente la función r = sqrt(R**2-z**2)
# Dado que el huevo tendrá simetría de revolución, "r" me representará el radio para cada punto a lo largo de z
# Como me interesa usar un espacio 3D mas grande que el huevo, para los valores de z fuera de la superficie del huevo: el radio=False 

R_xy = (x*x)+(y*y)
R_z = (R*R)-(z*z)
M = R_xy<=R_z

tf = time.time()

print("Tiempo: "+str(tf-t0))
print("Tamaño en MB: "+str(x.nbytes/1000000))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
scatter = ax.scatter3D(x, y, z, c=M, cmap='gray', s=0.001)
plt.show()


# Implementación del método de diferencias finitas para la ecuación de difusión de calor, usando una función

# Implementar la variación de la generación de calor en la yema

# Implementar el cambio de volumen de los componentes del huevo.


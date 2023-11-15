import numpy as np
import matplotlib.pylab as plt
import time

# Malla cúbica que contiene las condiciones iniciales del sistema, en el interior se encontrará el huevo.

# Generaré una matriz de nxnxn con los valores True y False que me van a definir los espacios: yema, albúmina, cáscara y exterior al huevo.

# Constantes de la simulación. 
								# En este caso elijo n=100
n = 100                                       			#500 es un buen número, consumo 250MB
m = 250								# m es el valor máximo y mínimo alcanzado por la matriz del espacio. En este caso 	
								# se define un espacio formado entre -25 y 25 en x,y,z. Unidades m x 10^-4 
d = np.int32(2*m/(n-1)) 					# Delta, distancia entre cada punto. d=5. Uso np.int32 para reducir procesamiento

# Constantes para la formación del huevo, en este caso el huevo es esférico. r=24mm
Rc = 240                                               		#Radio máximo de la cáscara                                      
Ra = 237 							#Radio máximo de la albúmina
Ry = 100							#Radio máximo de la yema
theta = np.radians(8)						#Este ángulo theta define la forma ovoide del huevo
g = np.tan(theta)

# Creación de 3 matrices que representan las coordenadas x,y,z para cada punto dentro de la matriz de espacio
t0 = time.time()

m1 = np.arange(-m,m+1,d,dtype=np.int32)                         # m1 es un arreglo de tamaño n que sirve de base para la construcción de las matrices x,y,z
x, y, z = np.meshgrid(m1,m1,m1, indexing='ij')                  # Creación de las matrices x,y,z. El orden de los vectores se toma de manera conveniente
      
####################################################################
# Yema - y
R_xy = (x*x)+(y*y)						#Distancia al eje z en x,y = 0
Ry_z = (Ry*Ry)-(z*z)						#Distancia a la yema para cada valor de z en x,y=0
My = R_xy<=Ry_z							#Región de la yema
####################################################################
# Albúmina - a
a_a = 242							#Constantes que definen la región de albúmina
b_a = 162

Ra_z = ((b_a+(g*z))*(b_a+(g*z)))*(1-((z*z)/(a_a*a_a)))		#Distancia a la albúmina para cada valor de z, en x,y=0
Ma = R_xy<=Ra_z							#Región interna a la superficie externa de la albúmina
####################################################################
# Cáscara - c
a_c = 245							#Constantes que definen la región cáscara
b_c = 165

Rc_z = ((b_c+(g*z))*(b_c+(g*z)))*(1-((z*z)/(a_c*a_c)))
Mc = R_xy<=Rc_z
#####################################################################
#Definición de las regiones:
Mc = Mc ^ Ma		#Región cáscara
Ma = Ma ^ My		#Región albúmina


#####################################################################
tf = time.time()
print("Tiempo: "+str(tf-t0))
print("Tamaño en MB: "+str(x.nbytes/1000000))

# Gráfica
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
scatter = ax.scatter3D(x, y, z, c=Ma*-1, cmap='gray', s=0.0005)                 #Gráfica de la cáscara
ax.view_init(0, 90)
plt.show()

##Prueba Git

# Implementación del método de diferencias finitas para la ecuación de difusión de calor, usando una función

# Implementar la variación de la generación de calor en la yema

# Implementar el cambio de volumen de los componentes del huevo.


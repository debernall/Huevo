#include <iostream>
#include <cmath>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;

int main () {
	ofstream outfile;									//Inicialización variable de documento
	outfile.open("Matriz.dat");								//Apertura de documento para escritura

	const double pi = M_PI;									//Pi
	int n = 100;										//Número de nodos en cada dirección de la matriz nxnxn
	double m = 25.0;									//Distancia máxima en cada eje desde el origen
	double dx = (2*m)/(n-1);								//Delta, separación entre cada nodo

	double Ry = 10.0;									//Radio máximo de la yema
	double g = tan(8*pi/180);								//Tangente en el punto medio del perfíl de huevo 
	double a_albu = 24.2;									//Longitud del semieje mayor albúmina
	double b_albu = 19.2;									//Longitud del semieje menor albúmina
	double a_casc = 24.6;									//Longitud del semieje mayor cáscara
	double b_casc = 19.6;									//Longitud del semieje menor cáscara

	vector<vector<vector<double>>> T(n, vector<vector<double>>(n, vector<double>(n)));	//Temperatura en t presente en cada nodo
	vector<vector<vector<double>>> Tf(n, vector<vector<double>>(n, vector<double>(n)));	//Temperatura en t futuro en cada nodo
	vector<vector<vector<double>>> Fo(n, vector<vector<double>>(n, vector<double>(n)));	//Número de Fourier para cada nodo

	double x,y,z;										//Inicialización de variables de longitud
	double r2,R2_casc,R2_albu,R2_yema;							//Inicialización de variables de radio
	
	double a = 0.000002;									//Difusividad térmica (PENDIENTE POR CAMBIAR - VALOR DE PRUEBA)
	double dt = 0.8*dx*dx/(6*a);								//Diferencial de tiempo calculado en base a la condición de estabilidad
	double t_muestreo = 1000000;								//Tiempo de muestreo máximo (PENDIENTE)
	int tf = int(t_muestreo/dt);								//Número máximo de intervalos de tiempo
	
	//CONDICIONES INICIALES
	double T0_Huevo=30;									//Temperatura inicial Huevo (TEMP DE PRUEBA)
	double T0_Aire=25;									//Temperatura inicial del Aire circundante
	
	for (int i=0;i<n;i++){
		for (int j=0; j<n;j++){
			for (int k=0; k<n; k++){
				x = (i*dx)-m;							//Posiciones en mm medidas desde el centro de la red de nodos
				y = (j*dx)-m;
				z = (k*dx)-m;
						
				r2 = x*x+z*z;								//Distancia de cada punto sobre sucesivos planos y al eje y
				R2_casc = (1-((y*y)/(a_casc*a_casc)))*(b_casc+(y*g))*(b_casc+(y*g));	//Distancia máxima de la cascara sobre sucesivos planos y al eje y
				R2_albu = (1-((y*y)/(a_albu*a_albu)))*(b_albu+(y*g))*(b_albu+(y*g));	//Distancia máxima de la albúmina sobre sucesivos planos y al eje y
				R2_yema = Ry*Ry - y*y;							//Distancia máxima de la yema sobre sucesivos planos y al eje y
				
				if (r2>R2_casc){						//Nodos externos al huevo
					T[i][j][k]=T0_Aire;						//Condición inicial de temperatura
					Fo[i][j][k]=a*dt/(dx*dx);					//Número de Fourier del nodo
				}
				else {
					if (r2<=R2_yema){					//Nodos pertenecientes a la yema
						T[i][j][k]=100;						//Condición inicial de temperatura (TEMPERATURA DE PRUEBA)
						Fo[i][j][k]=a*dt/(dx*dx);				//Número de Fourier del nodo
					}
					else if (r2<=R2_albu){					//Nodos pertenecientes a la albúmina
						T[i][j][k]=T0_Huevo;					//Condición inicial de temperatura
						Fo[i][j][k]=a*dt/(dx*dx);				//Número de Fourier del nodo
					}
					else {							//Nodos pertenecientes a la cáscara
						T[i][j][k]=T0_Huevo;					//Condición inicial de temperatura
						Fo[i][j][k]=a*dt/(dx*dx);				//Número de Fourier del nodo
					}	
				}
			}	
		}
	}	
	
	
	//EVOLUCIÓN TEMPORAL
	double Fo_i;								//Inicialización de variable Número de Fourier para reducir la extensión de la ecuación de calor
	double k = -0.0000001;
	double rmax_y = 19.0;
	double cont = 0;
	cout << Ry << "\n";
	for (int t=0; t<tf; t++){
		cont = cont +1;
		for (int i=0;i<n;i++){
			for (int j=0; j<n;j++){
				for (int k=0; k<n; k++){
					Fo_i = Fo[i][j][k];					//Almaceno el valor de Fo del nodo en una variable double para simplificar escritura
					x = (i*dx)-m;							//Posiciones en mm medidas desde el centro de la red de nodos
					y = (j*dx)-m;
					z = (k*dx)-m;
					
					
					
					//Ry = Ry * cbrt(1+(k*dt*(log(rmax_y)-log(4*pi*Ry*Ry*Ry/3))));
					
					r2 = x*x+z*z;								//Distancia de cada punto sobre sucesivos planos y al eje y
					R2_casc = (1-((y*y)/(a_casc*a_casc)))*(b_casc+(y*g))*(b_casc+(y*g));	//Distancia máxima de la cascara sobre sucesivos planos y al eje y
					R2_albu = (1-((y*y)/(a_albu*a_albu)))*(b_albu+(y*g))*(b_albu+(y*g));	//Distancia máxima de la albúmina sobre sucesivos planos y al eje y
					R2_yema = Ry*Ry - y*y;							//Distancia máxima de la yema sobre sucesivos planos y al eje y
					
					if (r2>R2_casc){						//Nodos externos al huevo
						Tf[i][j][k]=25;					//CONDICIÓN DE CONTORNO: El aire circundante se mantiene a temperatura constante
					}
					else {
						if (r2<=R2_yema){					//Nodos pertenecientes a la yema
							T[i][j][k]=100;
							Fo[i][j][k]=a*dt/(dx*dx);				//Número de Fourier del nodo
						}
						else if (r2<=R2_albu){					//Nodos pertenecientes a la albúmina
							Fo[i][j][k]=a*dt/(dx*dx);				//Número de Fourier del nodo
						}
						else {							//Nodos pertenecientes a la cáscara
							Fo[i][j][k]=a*dt/(dx*dx);				//Número de Fourier del nodo
						}
						//Temperatura futura por medio del método diferencias finitas - balance de energía
						//Tf[i][j][k]=Fo_i*((T[i][j][k]*((1/Fo_i)-6))+T[i+1][j][k]+T[i-1][j][k]+T[i][j+1][k]+T[i][j-1][k]+T[i][j][k+1]+T[i][j][k-1]);	
					}				
				}
			}
		}
		for (int i=0;i<n;i++){								//Almaceno la Temperatura futura calculada en Temperatura presente para el sig ciclo 
			for (int j=0; j<n;j++){
				for (int k=0; k<n; k++){
					T[i][j][k]=T[i][j][k];
				}
			}
		}
	}
	cout << dt << "\n";
	cout << Ry << "\n";
	cout << "pow " << cbrt(1+(k*dt*(log(rmax_y)-log(4*pi*Ry*Ry*Ry/3)))) << "\n";
	cout << "pow " << 1+(k*dt*(log(rmax_y)-log(4*pi*Ry*Ry*Ry/3))) << "\n";
	
	cout << Ry << "\n";
	cout << cont ;
	//IMPRESIÓN A ARCHIVO EXTERNO
	for (int i=0;i<n;i++){
			for (int j=0; j<n;j++){
				for (int k=0; k<n; k++){
					outfile << T[i][j][k] << "\n";				//Almaceno la última temperatura en un documento externo
				}
			}
		}


	outfile.close();									//Cierre del documento para escritura
	return 0;
}

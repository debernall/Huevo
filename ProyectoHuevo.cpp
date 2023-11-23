#include <iostream>
#include <cmath>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <algorithm>
#include <regex>

using namespace std;

int main () {

	//********************IMPORTE DE VARIABLES Y DATOS INICIALES****************************
	string linea;										//Variable string que me va a almacenar la info de las variables
  	fstream Datos ("DatosIniciales.dat");							//Importación de los datos iniciales
	int cont=0;
	map<string, double> variables;								//Creo un diccionario que va a almacenar las variables
												//Fuente: https://en.cppreference.com/w/cpp/container/map
	while ( getline (Datos,linea) ){							//Se recorre cada linea del documento
		regex exp_variable ("(.*)=");							//Expresión regular, busco lo que está despues del igual. () es agrupación,  					
		regex exp_valor ("=(.*)");							//. es cualquier caracter, * es un núm indeterminado de caracteres
												//Fuente: https://www.vishalchovatiya.com/regex-c/
		smatch variable;								//Aquí se va a almacenar el nombre de la variable
		smatch valor;									//Aquí se va a almacenar el valor de la variable
		regex_search(linea, variable, exp_variable);					//Se hace la busqueda en cada linea con la expresión regular
		regex_search(linea, valor, exp_valor);						//Fuente: https://en.cppreference.com/w/cpp/regex/regex_search
		variables[variable[1]]=stod(valor[1]);						//Almaceno cada variable en el diccionario de variables
    	}
	Datos.close();										//Cierre del archivo DatosIniciales

	//*************************DECLARACIÓN DE VARIABLES*************************************
	//Escalares del modelo
	const double pi = M_PI;									//Pi
	//Geometría del modelo
	int n = variables["n"];									//Número de nodos en cada dirección de la matriz nxnxn
	double m = variables["m"];								//Distancia máxima en cada eje desde el origen
	double dx = (2*m)/(n-1);								//Delta, separación entre cada nodo
	double Ry = variables["Ry"];								//Radio máximo de la yema
	double theta = variables["theta"];
	double g = tan(theta*pi/180);								//Tangente en el punto medio del perfíl de huevo 
	double a_albu = variables["a_albu"];							//Longitud del semieje mayor albúmina
	double b_albu = variables["b_albu"];							//Longitud del semieje menor albúmina
	double a_casc = variables["a_casc"];							//Longitud del semieje mayor cáscara
	double b_casc = variables["b_casc"];							//Longitud del semieje menor cáscara
	double x,y,z;										//Inicialización de variables de longitud
	double r2,R2_casc,R2_albu,R2_yema;							//Inicialización de variables de radio
	double rmax_y = variables["rmax_y"];							//Radio máximo de la yema
	//Propiedades térmicas	
	double a = variables["a"];								//Difusividad térmica (PENDIENTE POR CAMBIAR - VALOR DE PRUEBA)	
	double T0_Aire=variables["T0_Aire"];							//Temperatura inicial del Aire circundante
	double T0_Huevo=variables["T0_Huevo"];							//Temperatura inicial Huevo (TEMP DE PRUEBA)
	//Evolución temporal
	double Fo_i;										//Inicialización de variable Número de Fourier 
	double dt = 0.8*dx*dx/(6*a);								//Diferencial de tiempo calculado en base a la condición de estabilidad
	double t_muestreo = dt*variables["t"];							//Tiempo de muestreo máximo (PENDIENTE)
	int tf = int(t_muestreo/dt);								//Número máximo de intervalos de tiempo
	int t_imp[4] = {int(variables["t0"]),int(variables["t1"]),int(variables["t2"]),int(variables["t3"])};	//Tiempos que se desean graficar
	//Crecimiento del embrión
	double kk = variables["kk"];								//Tasa de crecimiento de la yema
	
	//Vectores del modelo
	vector<vector<vector< double >>> T(n, vector<vector< double >>(n, vector< double >(n)));//Temperatura en t presente en cada nodo
	vector<vector<vector<double>>> Tf(n, vector<vector<double>>(n, vector<double>(n)));	//Temperatura en t futuro en cada nodo
	vector<vector<vector<double>>> Fo(n, vector<vector<double>>(n, vector<double>(n)));	//Número de Fourier para cada nodo
		
	//*************************CONDICIONES INICIALES****************************************
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
	
	//****************************EVOLUCIÓN TEMPORAL****************************************
	for (int t=0; t<tf; t++){
		Ry = Ry * cbrt(1+(kk*dt*(log(rmax_y)-log(4*pi*Ry*Ry*Ry/3))));			//Modelo Gompertz del crecimiento de un conjunto de celulas, variable: Radio de Yema
		//*******************Cálculo de la temperatura en cada nodo********************
		for (int i=0;i<n;i++){
			for (int j=0; j<n;j++){
				for (int k=0; k<n; k++){
					Fo_i = Fo[i][j][k];					//Almaceno el valor de Fo del nodo en una variable double para simplificar escritura
					x = (i*dx)-m;						//Posiciones en mm medidas desde el centro de la red de nodos
					y = (j*dx)-m;
					z = (k*dx)-m;
					r2 = x*x+z*z;								//Distancia de cada punto sobre sucesivos planos y al eje y
					R2_casc = (1-((y*y)/(a_casc*a_casc)))*(b_casc+(y*g))*(b_casc+(y*g));	//Distancia máxima de la cascara sobre sucesivos planos y al eje y
					R2_albu = (1-((y*y)/(a_albu*a_albu)))*(b_albu+(y*g))*(b_albu+(y*g));	//Distancia máxima de la albúmina sobre sucesivos planos y al eje y
					R2_yema = Ry*Ry - y*y;							//Distancia máxima de la yema sobre sucesivos planos y al eje y

					if (r2>R2_casc){					//Nodos externos al huevo
						Tf[i][j][k]=25;					//CONDICIÓN DE CONTORNO: El aire circundante se mantiene a temperatura constante
					}
					else {
						if (r2<=R2_yema){				//Nodos pertenecientes a la yema
							T[i][j][k]=100;
							Fo[i][j][k]=a*dt/(dx*dx);		//Número de Fourier del nodo
						}
						else if (r2<=R2_albu){				//Nodos pertenecientes a la albúmina
							Fo[i][j][k]=a*dt/(dx*dx);		//Número de Fourier del nodo
						}
						else {						//Nodos pertenecientes a la cáscara
							Fo[i][j][k]=a*dt/(dx*dx);		//Número de Fourier del nodo
						}
						//Temperatura futura por medio del método diferencias finitas - balance de energía
						Tf[i][j][k]=Fo_i*((T[i][j][k]*((1/Fo_i)-6))+T[i+1][j][k]+T[i-1][j][k]+T[i][j+1][k]+T[i][j-1][k]+T[i][j][k+1]+T[i][j][k-1]);	
					}				
				}
			}
		}
		//******Se almacena la temperatura futura en la temperatura presente************	
		for (int i=0;i<n;i++){								 
			for (int j=0; j<n;j++){
				for (int k=0; k<n; k++){
					T[i][j][k]=Tf[i][j][k];
				}
			}
		}
		//*********Escritura a archivo externo para tiempos predeterminados**************
		for (int u=0; u<4;u++){	
			if (t==t_imp[u]){
				string Nombre = "Matriz";
				Nombre = Nombre + to_string(t_imp[u]) + ".dat";			//Me permite crear documentos independientes para cada tiempo solicitado
				ofstream outfile;						//Inicialización variable de documento
				outfile.open(Nombre);						//Apertura de documento para escritura				
				for (int i=0;i<n;i++){						
					for (int j=0; j<n;j++){
						for (int k=0; k<n; k++){
							outfile << T[i][j][k] << "\n";	
						}
					}
				}
				outfile.close();
			}
		}
	}									
	return 0;
}

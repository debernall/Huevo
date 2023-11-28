#include <iostream>
#include <cmath>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <algorithm>
#include <regex>
#include <ctime>

using namespace std;

//*********************************************PROYECTO HUEVO
vector<double> k_term(7,0);							//Conductividad térmica
vector<double> rho(7,0);							//Densidad
vector<double> cp(7,0);								//Capacidad calorífica
vector<double> a(7,0);								//Difusividad térmica
double dx;									
double dt;

//*********Función que calcula la temperatura para cada punto en cada instante de tiempo
double BalEnerg(double P0, double T0, double T1, double T2, double T3, double T4, double T5, double T6,	double Te1, double Te2, double Te3, double Te4, double Te5, double Te6, double Te7, double Te8, double P1, double P2, double P3, double P4, double P5, double P6, double P7, double P8){
	double Fo=k_term[P0]*dt/(rho[P0]*cp[P0]*dx*dx);				//Número de Fourier, me permite calcular faciltamente el valor adecuado de estabilidad
	double T;								//Variable para almacenar la temperatura calculada
	//bool parche=0;							//Variable booleana que me permite identificar si el nodo corresponde al parche de incubación
	//double K[8]={0,0,0,0,0,0,0,0};			
	//double P[8]={P1,P2,P3,P4,P5,P6,P7,P8};				
	
	//Estaba intentando implementar diferencias finitas aplicando solamente el térmico de conduccción para nodos que pertenecen a la frontera de cada región descrita pero no me funcionó porque requiere un dt demasiado pequeño y los resultados no se ajustaron a los esperados.
	
	//int VInt=0;
	//for (int e=0;e<8;e++){
	//	if (P[e]==P0){
	//		K[e]=k_term[P0];
	//		VInt=VInt+1;
	//	}
	//	else {
	//		K[e]=k_term[P[e]];
	//		if (P[e]==2 && P0==3){
	//			parche=1;
	//		}
	//	}
	//}
	if (P0<=1){								
		T=T0;
	}
//	else if(P0==6){
//		T=38;
//	}
	else{
		T = Fo*((T0*((1/Fo)-6))+T1+T2+T3+T4+T5+T6);
	}
	
/*	if (VInt>0 && P0>2){


		T=0;
		T=T+((K[0]+K[2]+K[4]+K[6])*(T1-T0));
		T=T+((K[1]+K[3]+K[5]+K[7])*(T2-T0));
		T=T+((K[0]+K[1]+K[4]+K[5])*(T3-T0));
		T=T+((K[2]+K[3]+K[6]+K[7])*(T4-T0));
		T=T+((K[0]+K[1]+K[2]+K[3])*(T5-T0));
		T=T+((K[4]+K[5]+K[6]+K[7])*(T6-T0));
		T = T0+((8*dt/(rho[P0]*cp[P0]*dx*VInt*4))*(T/dx));
		if (P0==6){
			T = T + (8*dt*368/(rho[P0]*cp[P0]*dx*VInt));
		}
		if (parche){
			T = T + ((8-VInt)*50000*dt/(rho[P0]*cp[P0]*dx*VInt));
		}
	}
	else {
		T=T0;
	}
*/		
	return T;
}

int main () {
	unsigned t0, t1;
	//********************IMPORTE DE VARIABLES Y DATOS INICIALES****************************
	string linea;										//Variable string que me va a almacenar la info de las variables
  	fstream Datos ("DatosInicialesCopia.dat");							//Importación de los datos iniciales
	map<string, double> variables;								//Creo un diccionario que va a almacenar las variables
												//Fuente: https://en.cppreference.com/w/cpp/container/map
	while ( getline (Datos,linea) ){							//Se recorre cada linea del documento
		regex exp_variable ("(.*)=");							//Expresión regular, busco lo que está despues del igual. () es agrupación,  					
		regex exp_valor ("=(.*)");							//. es cualquier caracter, * es un núm indeterminado de caracteres
												//Fuente: https://www.vishalchovatiya.com/regex-c/
		smatch variable;								//Aquí se va a almacenar el nombre de la variable
		smatch valor;									//Aquí se va a almacenar el valor de la variable
		if (regex_search(linea, variable, exp_variable)){				//Se hace la busqueda en cada linea con la expresión regular
			if(regex_search(linea, valor, exp_valor)){				//Fuente: https://en.cppreference.com/w/cpp/regex/regex_search
				variables[variable[1]]=stod(valor[1]);				//Almaceno cada variable en el diccionario de variables
    			}
    		}
    	}
	Datos.close();										//Cierre del archivo DatosIniciales

	//*************************DECLARACIÓN DE VARIABLES*************************************
	//Escalares del modelo
	const double pi = M_PI;									//Pi
	//Geometría del modelo
	int n = variables["n"];									//Número de nodos en cada dirección de la matriz nxnxn
	double m = variables["m"];								//Distancia máxima en cada eje desde el origen
	dx = (2*m)/(n-1);									//Delta, separación entre cada nodo
	double Ry = variables["Ry"];								//Radio máximo de la yema
	double theta = variables["theta"];
	double g = tan(theta*pi/180);								//Tangente en el punto medio del perfíl de huevo 
	double a_albu = variables["a_albu"];							//Longitud del semieje mayor albúmina
	double b_albu = variables["b_albu"];							//Longitud del semieje menor albúmina
	double a_casc = variables["a_casc"];							//Longitud del semieje mayor cáscara
	double b_casc = variables["b_casc"];							//Longitud del semieje menor cáscara
	double x,y,z;										//Inicialización de variables de longitud
	double r2,R2_casc,R2_albu,R2_yema,R_embr,R2_embr;					//Inicialización de variables de radio
	double Nodos[7]={0,0,0,0,0,0,0};							//En este arreglo me va a permitir contar los nodos de cada región 
	
	//Propiedades térmicas	
	k_term = {variables["k_aire"],variables["k_nido"],variables["k_parc"],variables["k_casc"],variables["k_albu"],variables["k_yema"],variables["k_embr"]};
	rho = {variables["rho_aire"],variables["rho_nido"],variables["rho_parc"],variables["rho_casc"],variables["rho_albu"],variables["rho_yema"],variables["rho_embr"]};
	cp = {variables["cp_aire"],variables["cp_nido"],variables["cp_parc"],variables["cp_casc"],variables["cp_albu"],variables["cp_yema"],variables["cp_embr"]};
	vector<double> a(7,0);									
	for (int i=0; i<7; i++){
		a[i]=k_term[i]/(rho[i]*cp[i]);	
	}
	double T0_Aire=variables["T0_Aire"];							//Temperatura inicial del Aire circundante
	double T0_Huevo=variables["T0_Huevo"];							//Temperatura inicial Huevo
	double T0_Parche=variables["T0_Parche"];						//Temperatura del parche de incubación
	double T0_Nido=variables["T0_Nido"];							//Temperatura del nido
	double T_Sup=variables["T0_Huevo"];							//Esta variable me ayuda a calcular el valor promedio de la temperatura del huevo
		
	//Evolución temporal
	dt = 0.95*dx*dx/(6*a[3]);								//Diferencial de tiempo calculado en base a la condición de estabilidad
	cout << dt << "\n";
	int tf = variables["t"];								//Número máximo de intervalos de tiempo
	int dt_imp = int(tf/50);								//Me permite imprimir en los archivos externos solo 500 valores de las variables de los tf periodos de tiempo
	int t_imp[4] = {0,33,66,100};								//Tiempos que se desean graficar 
	
	//Crecimiento del embrión
	double kk = variables["kk"];								//Tasa de crecimiento de la yema
	
	//Vectores del modelo
	vector<vector<vector< double >>> T(n, vector<vector< double >>(n, vector< double >(n)));//Temperatura en t presente en cada nodo
	vector<vector<vector<double>>> Tf(n, vector<vector<double>>(n, vector<double>(n)));	//Temperatura en t futuro en cada nodo
	vector<vector<vector<double>>> Fo(n, vector<vector<double>>(n, vector<double>(n)));	//Número de Fourier para cada nodo
	vector<double> TMuestreo(6,0);								//Almacena la temperatura media de cada región
	vector<double> VolMuestreo(6,0);							//Almacena el volumen de cada región
	vector<vector<vector< double >>> Pos(n, vector<vector< double >>(n, vector< double >(n)));	//Almacena un valor entre 0 y 6 que me define a que estructura corresponde el nodo
	
	//Archivo de muestreo de temperatura media y volumen
	ofstream outfileTemp;
	outfileTemp.open("TemperaturaMediaCopia.dat");
	outfileTemp << "Tiempo ; Temp_Aire ; Temp_Nido ; Temp_Parche; Temp_Cascara ; Temp_Albumina ; Temp_Yema ; Temp_Embrión\n";
	
	ofstream outfileVol;
	outfileVol.open("VolumenCopia.dat");
	outfileVol << "Tiempo ; Vol_Aire ; Vol_Nido ; Vol_Parche; Vol_Cascara ; Vol_Albumina ; Vol_Yema ; Vol_Embrión\n";
	
	
	//*************************CONDICIONES INICIALES****************************************
	for (int i=0;i<n;i++){
		
		for (int j=0; j<n;j++){
			for (int k=0; k<n; k++){
				x = (i*dx)-m;								//Posiciones en mm medidas desde el centro de la red de nodos
				y = (j*dx)-m;
				z = (k*dx)-m;
				r2 = x*x+z*z;								//Distancia de cada punto sobre sucesivos planos y al eje y
				R2_casc = (1-((y*y)/(a_casc*a_casc)))*(b_casc+(y*g))*(b_casc+(y*g));	//Distancia máxima de la cascara sobre sucesivos planos y al eje y
				R2_albu = (1-((y*y)/(a_albu*a_albu)))*(b_albu+(y*g))*(b_albu+(y*g));	//Distancia máxima de la albúmina sobre sucesivos planos y al eje y
				R2_yema = Ry*Ry - y*y;							//Distancia máxima de la yema sobre sucesivos planos y al eje y
				R_embr = 0.002;								//Radio inicial del embrión
				R2_embr = R_embr*R_embr;						//Radio^2 del embrión
				if (r2>R2_casc){							//Nodos externos al huevo				
					if (z>(b_casc*0.7)){						//Nodo parche incubación
						T[i][j][k]=T0_Parche;
						Pos[i][j][k]=2;
					}
					else if	(z<-(b_casc*0.8)){					//Nodo nido
						T[i][j][k]=T0_Nido;
						Pos[i][j][k]=1;
					}	
					else {								//Nodo aire
						T[i][j][k]=T0_Aire;
						Pos[i][j][k]=0;
					}
				}
				
				else {
					if (r2<=R2_embr){						//Nodos pertenecientes al embrión
						T[i][j][k]=T0_Huevo;					
						Pos[i][j][k]=6;				
					}
					else if (r2<=R2_yema){						//Nodos pertenecientes a la yema
						T[i][j][k]=T0_Huevo;					
						Pos[i][j][k]=5;
					}
					else if (r2<=R2_albu){						//Nodos pertenecientes a la albúmina
						T[i][j][k]=T0_Huevo;					
						Pos[i][j][k]=4;
					}
					else {								//Nodos pertenecientes a la cáscara
						T[i][j][k]=T0_Huevo;					
						Pos[i][j][k]=3;
					}	
				}
			}	
		}
	}	
	
	t0=clock();
	//****************************EVOLUCIÓN TEMPORAL****************************************
	bool impr;
	bool impr2;
	
	for (int t=0; t<tf; t++){
		
		if (t%dt_imp==0){impr=1;}else{impr=0;}							//Las variables impr e impr2 determinan si los datos de las variables de un periodo de tiempo dado se van a imprimir en archivo externo para graficar	
		if (t%(int(tf/3))==0){impr2=1;}else{impr2=0;}
		R_embr = R_embr * cbrt(1+(kk*dt*(log(4*pi*0.9*b_albu*b_albu*b_albu/3)-log(4*pi*R_embr*R_embr*R_embr/3))));										
		//Modelo Gompertz del crecimiento de un conjunto de celulas, variable: Radio del embrión
		//*******************Evolución espacial del huevo********************
		for (int i=0;i<n;i++){
			for (int j=0; j<n;j++){
				for (int k=0; k<n; k++){
					x = (i*dx)-m;								//Posiciones en mm medidas desde el centro de la red de nodos
					y = (j*dx)-m;
					z = (k*dx)-m;
					r2 = x*x+z*z;								//Distancia de cada punto sobre sucesivos planos y al eje y
					R2_casc = (1-((y*y)/(a_casc*a_casc)))*(b_casc+(y*g))*(b_casc+(y*g));	//Distancia máxima de la cascara sobre sucesivos planos y al eje y
					R2_albu = (1-((y*y)/(a_albu*a_albu)))*(b_albu+(y*g))*(b_albu+(y*g));	//Distancia máxima de la albúmina sobre sucesivos planos y al eje y
					R2_embr = (1-((y*y)/(R_embr*R_embr/(0.793388*0.793388))))*(R_embr+(y*g))*(R_embr+(y*g));
					R2_yema = Ry*Ry - y*y;							//Distancia máxima de la yema sobre sucesivos planos y al eje y
					if (r2>R2_casc){											
						if (z>(b_casc*0.7)){
							
								T[i][j][k]=T0_Aire;
								Pos[i][j][k]=0;
								Nodos[0]=Nodos[0]+1;
							
						}
						else if	(z<-(b_casc*0.8)){
							Pos[i][j][k]=1;
							Nodos[1]=Nodos[1]+1;
							}
						else {
							Pos[i][j][k]=0;
							Nodos[0]=Nodos[0]+1;
							}
					}
					else {		
						if (r2<=R2_embr){				
							Pos[i][j][k]=6;
							Nodos[6]=Nodos[6]+1;
						}						
						else if (r2<=R2_yema){			
							Pos[i][j][k]=5;
							Nodos[5]=Nodos[5]+1;
						}
						else if (r2<=R2_albu){			
							Pos[i][j][k]=4;
							Nodos[4]=Nodos[4]+1;
						}
						else {					
							Pos[i][j][k]=3;
							Nodos[3]=Nodos[3]+1;
						}	
					}
				}	
			}
		}

		//*******************Evolución de la temperatura del huevo********************
		for (int i=0;i<n;i++){
			for (int j=0; j<n;j++){
				for (int k=0; k<n; k++){		
					if (i>0 && i<n-1 && j>0 && j<n-1 && k>0 && k<n-1){			//Selecciona nodos que internos a la frontera
						Tf[i][j][k]=BalEnerg(Pos[i][j][k],T[i][j][k],
								T[i+1][j][k],T[i-1][j][k],
								T[i][j+1][k],T[i][j-1][k],
								T[i][j][k+1],T[i][j][k-1],
								T[i+1][j+1][k+1],T[i-1][j+1][k+1],
								T[i+1][j-1][k+1],T[i-1][j-1][k+1],
								T[i+1][j+1][k-1],T[i-1][j+1][k-1],
								T[i+1][j-1][k-1],T[i-1][j-1][k-1],
								Pos[i+1][j+1][k+1],Pos[i-1][j+1][k+1],
								Pos[i+1][j-1][k+1],Pos[i-1][j-1][k+1],
								Pos[i+1][j+1][k-1],Pos[i-1][j+1][k-1],
								Pos[i+1][j-1][k-1],Pos[i-1][j-1][k-1]);
					}
					else if(Pos[i][j][k]==0){						//Condiciones de frontera abierta
						Tf[i][j][k]=T0_Aire;
					}
					else if(Pos[i][j][k]==1){
						Tf[i][j][k]=T0_Nido;
					}
					else if(Pos[i][j][k]==2){
						Tf[i][j][k]=T0_Parche;
					}		
					TMuestreo[Pos[i][j][k]]=TMuestreo[Pos[i][j][k]]+Tf[i][j][k];		//Se almacena la temperatura de cada nodo para calcular la temperatura promedio de cada estructura
				}
			}
		}
		T_Sup=((TMuestreo[3]/Nodos[3])+(TMuestreo[4]/Nodos[4])+(TMuestreo[6]/Nodos[6]))/3;		//Temperatura promedio del huevo, no tengo en cuenta la estructura yema ya que con el tiempo esta estructura desaparece y por tanto no se tendría valor para la temperatura de esta estructura

		//******Se almacena la temperatura futura en la temperatura presente************	
		for (int i=0;i<n;i++){								 
			for (int j=0; j<n;j++){
				for (int k=0; k<n; k++){
					T[i][j][k]=Tf[i][j][k];
				}
			}
		}
		
		//*********Escritura a archivo externo para tiempos predeterminados**************
		if (impr2){
			string Nombre = "Matriz";
			Nombre = Nombre + to_string(t_imp[int(t/(int(tf/3)))]) + "Copia.dat";//Me permite crear documentos independientes para cada tiempo solicitado
			ofstream outfile;						//Inicialización variable de documento
			outfile.open(Nombre);						//Apertura de documento para escritura				
			for (int i=0;i<n;i++){						
				for (int j=0; j<n;j++){
					for (int k=0; k<n; k++){
						if (Pos[i][j][k]>2){
							outfile << T[i][j][k] << "\n";
						}
						else {
							outfile << 0 << "\n";
						}					
					}
				}
			}
			outfile.close();
		}		
		
		//*************Temperatura media y volumen de cada región***********************************	
		if (impr){
			outfileTemp << t*dt/60 << " ; ";
			outfileVol << t*dt/60 << " ; ";
			for (int w=0;w<7;w++){
				TMuestreo[w]=TMuestreo[w]/Nodos[w];
				VolMuestreo[w]=Nodos[w]*dx*dx*dx;
				outfileTemp << TMuestreo[w];
				outfileVol << VolMuestreo[w];
				if (w<6){
					outfileTemp << " ; ";
					outfileVol << " ; ";
				}
				else {
					outfileTemp << "\n";
					outfileVol << "\n";
				}
				VolMuestreo[w]=0;								
			}
		}
		for (int w=0;w<7;w++){
			Nodos[w]=0;
			TMuestreo[w]=0;		
		} 		
  	}	
	outfileTemp.close();
	outfileVol.close();
	t1=clock();
	int time=int((double(t1-t0)/CLOCKS_PER_SEC));
	cout << time << "\n";
							
	return 0;
}

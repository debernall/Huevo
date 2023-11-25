#include <iostream>
#include <cmath>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <algorithm>
#include <regex>
#include <ctime>

using namespace std;

double BalEnerg(double P0, double T0, double T1, double T2, double T3, double T4, double T5, double T6, double P1, double P2, double P3, double P4, double P5, double P6){
	double P[6] = {P1,P2,P3,P4,P5,P6};
	int Nod_int=0;
	for (int i=0;i<6;i++){
		if (P0==P[i]){
			Nod_int = Nod_int +1;
		}
	}
	//Fo*((T0*((1/Fo)-6))+T1+T2+T3+T4+T5+T6);
	if (P0==6){T0=100;}
	
	return T0;
}

int main () {
	unsigned t0, t1;
	
	//********************IMPORTE DE VARIABLES Y DATOS INICIALES****************************
	string linea;										//Variable string que me va a almacenar la info de las variables
  	fstream Datos ("DatosIniciales.dat");							//Importación de los datos iniciales
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
	double dx = (2*m)/(n-1);								//Delta, separación entre cada nodo
	double Ry = variables["Ry"];								//Radio máximo de la yema
	double theta = variables["theta"];
	double g = tan(theta*pi/180);								//Tangente en el punto medio del perfíl de huevo 
	double a_albu = variables["a_albu"];							//Longitud del semieje mayor albúmina
	double b_albu = variables["b_albu"];							//Longitud del semieje menor albúmina
	double a_casc = variables["a_casc"];							//Longitud del semieje mayor cáscara
	double b_casc = variables["b_casc"];							//Longitud del semieje menor cáscara
	double x,y,z;										//Inicialización de variables de longitud
	double r2,R2_casc,R2_albu,R2_yema,R_embr,R2_embr;							//Inicialización de variables de radio
	double rmax_y = variables["rmax_y"];							//Radio máximo de la yema
	//Propiedades térmicas	
	double k = variables["k"];
	double rho = variables["rho"];
	double cp = variables["cp"];
	double a = k/(rho*cp);								//Difusividad térmica (PENDIENTE POR CAMBIAR - VALOR DE PRUEBA)	
	double T0_Aire=variables["T0_Aire"];							//Temperatura inicial del Aire circundante
	double T0_Huevo=variables["T0_Huevo"];							//Temperatura inicial Huevo (TEMP DE PRUEBA)
	double Nodos[7]={0,0,0,0,0,0,0};
	//Evolución temporal
	double Fo_i;										//Inicialización de variable Número de Fourier 
	double dt = 0.95*dx*dx/(6*a);								//Diferencial de tiempo calculado en base a la condición de estabilidad
	
	double t_m = dt*variables["t"];							//Tiempo de muestreo máximo (PENDIENTE)
	int tf = variables["t"];								//Número máximo de intervalos de tiempo
	
	//int t_muestreo[100]={0};
	//for (int i=0;i<100;i++){t_muestreo[i]=int(i*tf/100);}
	int dt_imp = int(tf/100);
	int t_imp[4] = {0,33,66,100};								//Tiempos que se desean graficar
//	int t_imp2[4] = {0,int(t_imp[1]*tf/100),int(t_imp[2]*tf/100),tf-1};
	
	//Crecimiento del embrión
	double kk = variables["kk"];								//Tasa de crecimiento de la yema
	double a_embr;
	double b_embr;
	double VolHuevo;
	int t_parche=int(5*60/dt);
	double f_parche=int(0.6*t_parche);
	
	//Vectores del modelo
	vector<vector<vector< double >>> T(n, vector<vector< double >>(n, vector< double >(n)));//Temperatura en t presente en cada nodo
	vector<vector<vector<double>>> Tf(n, vector<vector<double>>(n, vector<double>(n)));	//Temperatura en t futuro en cada nodo
	vector<vector<vector<double>>> Fo(n, vector<vector<double>>(n, vector<double>(n)));	//Número de Fourier para cada nodo
	vector<double> TMuestreo(6,0);				//Almacena la temperatura media de cada región
	vector<double> VolMuestreo(6,0);				//Almacena la temperatura media de cada región
	vector<vector<vector< double >>> Pos(n, vector<vector< double >>(n, vector< double >(n)));
	
	//Archivo de muestreo de temperatura media
	ofstream outfileTemp;
	outfileTemp.open("TemperaturaMedia.dat");
	outfileTemp << "Tiempo ; Temp_Aire ; Temp_Nido ; Temp_Parche; Temp_Cascara ; Temp_Albumina ; Temp_Yema ; Temp_Embrión\n";
	
	ofstream outfileVol;
	outfileVol.open("Volumen.dat");
	outfileVol << "Tiempo ; Vol_Aire ; Vol_Nido ; Vol_Parche; Vol_Cascara ; Vol_Albumina ; Vol_Yema ; Vol_Embrión\n";
	
	//*************************CONDICIONES INICIALES****************************************
	
	for (int i=0;i<n;i++){
		
		for (int j=0; j<n;j++){
			for (int k=0; k<n; k++){
				x = (i*dx)-m;						//Posiciones en mm medidas desde el centro de la red de nodos
				y = (j*dx)-m;
				z = (k*dx)-m;
				r2 = x*x+z*z;								//Distancia de cada punto sobre sucesivos planos y al eje y
				R2_casc = (1-((y*y)/(a_casc*a_casc)))*(b_casc+(y*g))*(b_casc+(y*g));	//Distancia máxima de la cascara sobre sucesivos planos y al eje y
				R2_albu = (1-((y*y)/(a_albu*a_albu)))*(b_albu+(y*g))*(b_albu+(y*g));	//Distancia máxima de la albúmina sobre sucesivos planos y al eje y
				R2_yema = Ry*Ry - y*y;						//Distancia máxima de la yema sobre sucesivos planos y al eje y
				R_embr = 0.001;
				R2_embr = R_embr*R_embr;
				if (r2>R2_casc){					//Nodos externos al huevo				
					if (z>(b_casc*0.7)){
						T[i][j][k]=40;
						Pos[i][j][k]=2;
					}
					else if	(z<-(b_casc*0.8)){
						T[i][j][k]=29;
						Pos[i][j][k]=1;
					}
					else {
						T[i][j][k]=T0_Aire;
						Pos[i][j][k]=0;
					}
				}
				
				else {
					if (r2<=R2_embr){					//Nodos pertenecientes a la yema
						T[i][j][k]=T0_Huevo;					//Condición inicial de temperatura (TEMPERATURA DE PRUEBA)
						Pos[i][j][k]=6;				
					}
					else if (r2<=R2_yema){					//Nodos pertenecientes a la yema
						T[i][j][k]=T0_Huevo;					//Condición inicial de temperatura (TEMPERATURA DE PRUEBA)
						Pos[i][j][k]=5;
					}
					else if (r2<=R2_albu){					//Nodos pertenecientes a la albúmina
						T[i][j][k]=T0_Huevo;					//Condición inicial de temperatura
						Pos[i][j][k]=4;
					}
					else {							//Nodos pertenecientes a la cáscara
						T[i][j][k]=T0_Huevo;					//Condición inicial de temperatura
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
		if (t%dt_imp==0){impr=1;}else{impr=0;}
		if (t%(int(tf/3))==0){impr2=1;}else{impr2=0;}
		R_embr = R_embr * cbrt(1+(kk*dt*(log(4*pi*rmax_y*rmax_y*rmax_y/3)-log(4*pi*R_embr*R_embr*R_embr/3))));										//Modelo Gompertz del crecimiento de un conjunto de celulas, variable: Radio de Yema
		
		//*******************Evolución espacial del huevo********************
		for (int i=0;i<n;i++){
			for (int j=0; j<n;j++){
				for (int k=0; k<n; k++){
					x = (i*dx)-m;						//Posiciones en mm medidas desde el centro de la red de nodos
					y = (j*dx)-m;
					z = (k*dx)-m;
					r2 = x*x+z*z;								//Distancia de cada punto sobre sucesivos planos y al eje y
					R2_casc = (1-((y*y)/(a_casc*a_casc)))*(b_casc+(y*g))*(b_casc+(y*g));	//Distancia máxima de la cascara sobre sucesivos planos y al eje y
					R2_albu = (1-((y*y)/(a_albu*a_albu)))*(b_albu+(y*g))*(b_albu+(y*g));	//Distancia máxima de la albúmina sobre sucesivos planos y al eje y
					R2_embr = (1-((y*y)/(R_embr*R_embr/(0.793388*0.793388))))*(R_embr+(y*g))*(R_embr+(y*g));
					R2_yema = Ry*Ry - y*y;							//Distancia máxima de la yema sobre sucesivos planos y al eje y
					
					if (r2>R2_casc){					//Nodos externos al huevo				
						if (z>(b_casc*0.7)){
							if ((t%t_parche)<f_parche){
								T[i][j][k]=40;
								Pos[i][j][k]=2;	
								Nodos[2]=Nodos[2]+1;	
								}
							else{
								T[i][j][k]=T0_Aire;
								Pos[i][j][k]=0;
								Nodos[0]=Nodos[0]+1;
							}
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
						if (r2<=R2_embr){				//Nodos pertenecientes a la yema
							Pos[i][j][k]=6;
							Nodos[6]=Nodos[6]+1;
						}						
						else if (r2<=R2_yema){				//Nodos pertenecientes a la yema
							Pos[i][j][k]=5;
							Nodos[5]=Nodos[5]+1;
						}
						else if (r2<=R2_albu){				//Nodos pertenecientes a la albúmina
							Pos[i][j][k]=4;
							Nodos[4]=Nodos[4]+1;
						}
						else {						//Nodos pertenecientes a la cáscara
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
					if (i>0 && i<n-1 && j>0 && j<n-1 && k>0 && k<n-1){
						Tf[i][j][k]=BalEnerg(Pos[i][j][k],T[i][j][k],
								T[i+1][j][k],T[i-1][j][k],
								T[i][j+1][k],T[i][j-1][k],
								T[i][j][k+1],T[i][j][k-1],
								Pos[i+1][j][k],Pos[i-1][j][k],
								Pos[i][j+1][k],Pos[i][j-1][k],
								Pos[i][j][k+1],Pos[i][j][k-1]);
					}
					else {Tf[i][j][k]=T0_Aire;}		
					TMuestreo[Pos[i][j][k]]=TMuestreo[Pos[i][j][k]]+Tf[i][j][k];
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
		if (impr2){
			string Nombre = "Matriz";
			Nombre = Nombre + to_string(t_imp[int(t/(int(tf/3)))]) + ".dat";			//Me permite crear documentos independientes para cada tiempo solicitado
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
		
		//*************Temperatura media de cada región***********************************	
		if (impr){
			outfileTemp << t*dt/3600 << " ; ";
			outfileVol << t*dt/3600 << " ; ";
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

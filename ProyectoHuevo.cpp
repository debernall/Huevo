#include <iostream>
#include <cmath>
#include <stdio.h>
#include <fstream>

using namespace std;

int main () {
	ofstream outfile;
	outfile.open("Matriz.dat");

	const double pi = M_PI;					//Pi
	int n = 50;						//Tamaño n^3
	double m = 25.0;					//Distancia máxima en cada eje desde el origen
	double d = (2*m)/(n-1);					//Delta, ancho de cada cuadrícula

	double Rc = 24.0;					//Radio máximo de la cáscara
	double Ra = 23.7;					//Radio máximo de la albúmina
	double Ry = 10.0;					//Radio máximo de la yema

	double g = tan(8*pi/180);				//Tangente en el punto medio del perfíl de huevo 

	double A[n][n][n];	

	double x;
	double y;
	double z;
	double r2;
	double R2_c;
	double R2_a;
	double R2_y;
	double a_a = 24.2;					//Longitud del semieje mayor albúmina
	double b_a = 16.2;					//Longitud del semieje menor albúmina
	double a_c = 24.5;					//Longitud del semieje mayor cáscara
	double b_c = 16.5;					//Longitud del semieje menor cáscara
	double alpha_y = 1;					//Difusividad términa de la yema
	double alpha_a = 2;					//Difusividad térmica de la albúmina
	double alpha_c = 3;					//Difusividad térmica de la cáscara

	for (int i=0;i<n;i++){
		for (int j=0; j<n;j++){
			for (int k=0; k<n; k++){
				x = (i*d)-m;
				y = (j*d)-m;
				z = (k*d)-m;
						
				r2 = x*x+y*y;
				R2_c = (1-((z*z)/(a_c*a_c)))*(b_c+(z*g))*(b_c+(z*g));
				R2_a = (1-((z*z)/(a_a*a_a)))*(b_a+(z*g))*(b_a+(z*g));
				R2_y = Ry*Ry - z*z;
				if (r2>R2_c){
					A[i][j][k]=0;
				}
				else {
					if (r2<=R2_y){
						A[i][j][k]=0;
					}
					else if (r2<R2_a){
						A[i][j][k]=200;	
					}
					else {
						A[i][j][k]=300;
					}	
				}
				outfile << A[i][j][k] << "\n";
			}	
		}
	}
	outfile.close();
	return 0;
}

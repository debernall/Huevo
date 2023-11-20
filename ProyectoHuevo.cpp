#include <iostream>
#include <cmath>
#include <stdio.h>
#include <fstream>

using namespace std;

int main () {
	ofstream outfile;
	outfile.open("Matriz.dat");

	const double pi = M_PI;					//Pi
	int n = 110;						//Tamaño n^3
	double m = 25.0;					//Distancia máxima en cada eje desde el origen
	double d = (2*m)/(n-1);					//Delta, ancho de cada cuadrícula

	double Ry = 10.0;					//Radio máximo de la yema

	double g = tan(8*pi/180);				//Tangente en el punto medio del perfíl de huevo 

	int A[n][n][n];	

	double x;
	double y;
	double z;
	double r2;
	double R2_c;
	double R2_a;
	double R2_y;
	double a_a = 24.2;					//Longitud del semieje mayor albúmina
	double b_a = 16.2;					//Longitud del semieje menor albúmina
	double a_c = 24.6;					//Longitud del semieje mayor cáscara
	double b_c = 16.6;					//Longitud del semieje menor cáscara
	double alpha_y = 1;					//Difusividad términa de la yema
	double alpha_a = 2;					//Difusividad térmica de la albúmina
	double alpha_c = 3;					//Difusividad térmica de la cáscara
	int error=0;
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
					else if (r2<=R2_a){
						A[i][j][k]=0;	
					}
					else if (r2<=R2_c){
						A[i][j][k]=1;
					}	
					else {
						error = error +1;				
					}
				}
				outfile << A[i][j][k] << "\n";
			}	
		}
	}
	
	cout << " error" << error << "\n";
	
	
	int ss[n][n][n];
	int s1=0;
	int s2=0;
	int s3=0;
	int s4=0;
	int s5=0;
	int s6=0;
	int s0=0;
	int se=0;
	int s1n=0;
	int s2n=0;
	int s3n=0;
	int s4n=0;
	int s5n=0;
	int s6n=0;
	int s0n=0;
	int sen=0;
	int cont = 0;
	int ttt = 0;
	
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
					ttt =0;
				}
				else {
					if (r2<=R2_y){
						ttt = 0;
					}
					else if (r2<R2_a){
						ttt =0;	
						
					}
					else {
						ttt =0;
						cont = cont +1;
						ss[i][j][k] = A[i+1][j][k]+A[i-1][j][k]+A[i][j+1][k]+A[i][j-1][k]+A[i][j][k+1]+A[i][j][k-1];
						if (ss[i][j][k]==0){s0=s0+1;}
						else if (ss[i][j][k]==1){s1=s1+1;}
						else if (ss[i][j][k]==2){s2=s2+1;}
						else if (ss[i][j][k]==3){s3=s3+1;}
						else if (ss[i][j][k]==4){s4=s4+1;}
						else if (ss[i][j][k]==5){s5=s5+1;}
						else if (ss[i][j][k]==6){s6=s6+1;}
						else if (ss[i][j][k]==-1){s1n=s1n+1;}
						else if (ss[i][j][k]==-2){s2n=s2n+1;}
						else if (ss[i][j][k]==-3){s3n=s3n+1;}
						else if (ss[i][j][k]==-4){s4n=s4n+1;}
						else if (ss[i][j][k]==-5){s5n=s5n+1;}
						else if (ss[i][j][k]==-6){s6n=s6n+1;}
						else {se=se+1;}
					}	
				}
			}	
		}
	}
	
	/*
	for (int i=1;i<n-1;i++){
		for (int j=1; j<n-1;j++){
			for (int k=1; k<n-1; k++){
				cont = cont +1;
				ss[i][j][k] = (A[i+1][j][k]+A[i-1][j][k]+A[i][j+1][k]+A[i][j-1][k]+A[i][j][k+1]+A[i][j][k-1])*A[i][j][k];
				if (ss[i][j][k]==0){s0=s0+1;}
				else if (ss[i][j][k]==1){s1=s1+1;}
				else if (ss[i][j][k]==2){s2=s2+1;}
				else if (ss[i][j][k]==3){s3=s3+1;}
				else if (ss[i][j][k]==4){s4=s4+1;}
				else if (ss[i][j][k]==5){s5=s5+1;}
				else if (ss[i][j][k]==6){s6=s6+1;}
				else if (ss[i][j][k]==-1){s1n=s1n+1;}
				else if (ss[i][j][k]==-2){s2n=s2n+1;}
				else if (ss[i][j][k]==-3){s3n=s3n+1;}
				else if (ss[i][j][k]==-4){s4n=s4n+1;}
				else if (ss[i][j][k]==-5){s5n=s5n+1;}
				else if (ss[i][j][k]==-6){s6n=s6n+1;}
				else {se=se+1;}
			}
		}
	}
	*/
	cout << " s0=" << s0 << "\n s1=" << s1 << "\n s2=" << s2 << "\n s3=" << s3 << "\n s4=" << s4 << "\n s5=" << s5 << "\n s6=" << s6 << "\n se=" << se; 
	cout << "\n s1n=" << s1n << "\n s2n=" << s2n << "\n s3n=" << s3n << "\n s4n=" << s4n << "\n s5n=" << s5n << "\n s6n=" << s6n ; 
	cout << "\n cont=" << cont;
	outfile.close();
	return 0;
}

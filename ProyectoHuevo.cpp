#include <iostream>
#include <cmath>
#include <stdio.h>
#include <fstream>

using namespace std;

int main () {
	ofstream outfile;
	outfile.open("Matriz.dat");

	const double pi = M_PI;
	int n = 100;
	double m = 25.0;
	double d = (2*m)/(n-1);

	double Rc = 24.0;
	double Ra = 23.7;
	double Ry = 10.0;

	double g = tan(8*pi/180);

	double A[n][n][n];

	double x;
	double y;
	double z;
	double r2;
	double R2;
	double a = 24.5;
	double b = 16.5;

	for (int i=0;i<n;i++){
		for (int j=0; j<n;j++){
			for (int k=0; k<n; k++){
				x = (i*d)-m;
				y = (j*d)-m;
				z = (k*d)-m;
						
				r2 = x*x+y*y;
				R2 = (1-((z*z)/(a*a)))*(b+(z*g))*(b+(z*g));

				if (r2>R2){
					A[i][j][k]=0;
				}
				else {
					A[i][j][k]=1;	
				}
				outfile << A[i][j][k] << "\n";
			}	
		}
	}
	outfile.close();
	return 0;
}

#include <iostream>
#include <cmath>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;

int main () {
	string line;
  	ifstream Datos ("DatosIniciales.dat");
	if (Datos.is_open())
	  {
	    while ( getline (Datos,line) )
	    {
	      cout << line << '\n';
	    }
	    Datos.close();
	  }

	  else cout << "Unable to open file"; 
		
	return 0;
}

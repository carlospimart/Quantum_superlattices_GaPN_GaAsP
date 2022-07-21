#include <iostream>
#include <vector>
#include <cmath>
using namespace std;
#define NL '\n'

//*********************************************************************************
//Funcion para el Calculo del espesor relativo de una superred en funcion de x e y 
//*********************************************************************************
double REL_Func(float a_binario1_a, float a_binario1_b, float a_binario2_a, float a_binario2_b, double x, double y){
                 
				    double REL;
				    float aSi=5.431;
				    REL=((aSi-(a_binario2_a*(1-y)+y*a_binario2_b)))/((a_binario1_a*x+(1-x)*a_binario1_b)-(a_binario2_a*(1-y)+y*a_binario2_b));
					   
					return REL;
				

}

//***********************************************************************************
//Programa para el Calculo del gap mediante la ley de vegard para un solo valor de x
//**********************************************************************************
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;
#define NL '\n'

double vegard_law_func(float EgrBinary_1, float EgrBinary_2, float CTernary_r, double x){
				  
				    double EgrTernary_c;
					EgrTernary_c=EgrBinary_1*x+EgrBinary_2*(1-x)-CTernary_r*x*(1-x);
					
					return EgrTernary_c;
				  
				
}

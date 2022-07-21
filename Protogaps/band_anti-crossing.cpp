#include <iostream>
#include <vector>
#include <cmath>
using namespace std;
#define NL '\n'

//*****************************************************************************
//Funcion para el Calculo del gap mediante BAC model para un solo valor de y:
//*******************************************************************************
double BAC_Model_func(float Em, float En, float Cmn, double y){
                 
				    double Egn_1;
				    
					Egn_1=0.5*((Em+En)-pow((pow((Em-En),2)+4*(y)*pow(Cmn,2)),0.5));
					   
					return Egn_1;
				

}

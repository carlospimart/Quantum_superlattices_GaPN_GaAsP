//************************************************************************************************************************************
//Programa para calcular la diferencia de la banda de conducci�n mas baja menos la banda de valencia mas alta en una superrred tipo II
//************************************************************************************************************************************
#include <iostream>
#include <vector>
#include "Calculator_Band_offset_General.cpp"
#include <cmath>
using namespace std;
#define NL '\n'
#include <fstream> 

int main(){
float Cmn_GaP=3.05;  //todo en eV
float Em_GaPN=2.78; //gap en el centro de la zona de Brillouin. El gap pasa de indirecto (gamma-X) a directo (gamma) cuando se añade nitrógeno.
float En_GaPN=2.18;   //referencias
double x;
double y;
double y_for_rel;
double protogap;


//******************************************************///
//**************parametros de red:**********************///
//******************************************************///
float aSi=5.431; //Amstrongs - Lattice parameter for silicon
float aGaAs=5.65325; //Amstrongs - lattice parametre a for GaAs
float aGaP=5.4505; //Amstrongs - lattice parametre a for GaP
float aGaN=4.52;//Amstrongs -  lattice parametre a for GaN
float aInP=5.8687;//Amstrongs -  lattice parametre a for InP. Ioffe
float aInN=5.004;//Amstrongs -  lattice parametre a for InN
float aGaIn=5.004;//Amstrongs -  lattice parametre a for GaIn


ofstream f_out("y_x_protogap.dat", ios::out);// archivo que guarda un mapa de protogaps para un valor de x e y dados	
ofstream f_out_1("rel_thick_for_x_y.dat", ios::out);//archivo declarado para guardar una recta de x e y para un espesor relativo constante
	
float rel_thick=0.8;

	
	for(x=0;x<0.61;x+=0.01){
	    
		y_for_rel=((aSi-rel_thick*(aGaAs*x+(1-x)*aGaP))/((1-rel_thick)*(aGaN-aGaP))) - (aGaP/(aGaN-aGaP));
	    
		f_out_1<<x<<", "<<y_for_rel<<NL;
		
	 for(y=0 ;y<0.151;y+=0.001){
	 	
          protogap  = BAC_Model_func(Em_GaPN, En_GaPN, Cmn_GaP, y) - Calculator_Band_offset_General(x, y, "VBO", "GaAsP/GaPN");
	      
		  cout<<"x: "<<x<<", y: "<<y<<NL;
		  cout<<"E_GAPN: "<<BAC_Model_func(Em_GaPN, En_GaPN, Cmn_GaP, y)<<", VBO: "<<Calculator_Band_offset_General(x, y, "VBO", "GaAsP/GaPN")<<", Protogap: "<<protogap<<NL;
		  
		  f_out<<x<<", "<<y<<", "<<protogap<<NL;
		  
		  
       }
	}
	
}

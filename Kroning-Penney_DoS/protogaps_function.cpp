//************************************************************************************************************************************
//Programa para calcular la diferencia de la banda de conducci�n mas baja menos la banda de valencia mas alta en una superrred tipo II
//************************************************************************************************************************************
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;
#define NL '\n'
#include <fstream> 

double protogaps_function(double x, double y, char *SL_type){
float Cmn_GaP=3.05;  //todo en eV
float Em_GaPN=2.78; //gap en el centro de la zona de Brillouin. El gap pasa de indirecto (gamma-X) a directo (gamma) cuando se añade nitrógeno.
float En_GaPN=2.18;   //referencias

double protogap;
double AntiProtogap;

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



	 	
	      protogap  = BAC_Model_func(Em_GaPN, En_GaPN, Cmn_GaP, y) - Calculator_Band_offset_General(x, y, "VBO", SL_type);
	      
	      AntiProtogap = BAC_Model_func(Em_GaPN, En_GaPN, Cmn_GaP, y) + Calculator_Band_offset_General(x, y, "CBO", SL_type);
	      
		  //protogap  = BAC_Model_func(Em_GaPN, En_GaPN, Cmn_GaP, y) - Band_Offset_GaInP_GaPN(x, y, "VBO");
		  
		  
		  cout<<"E_GAPN: "<<BAC_Model_func(Em_GaPN, En_GaPN, Cmn_GaP, y)<<", Protogap: "<<protogap<<NL;
		  
		  cout<<"Anti-Protogap: "<<AntiProtogap<<NL;
		  
return protogap;
	
}

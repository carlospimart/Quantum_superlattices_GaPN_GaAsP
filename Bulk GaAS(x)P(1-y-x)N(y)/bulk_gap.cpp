#include <iostream>
#include <vector>
#include <cmath>
using namespace std;
#define NL '\n'
#include <fstream> 


//*****************************************************************************
//Función para el Cálculo del gap mediante BAC model para un solo valor de y:
//*******************************************************************************
double BAC_Model(float Em, float En, float Cmn, double y){
                 
				    double Egn_1;
					Egn_1=0.5*((Em+En)-pow((pow((Em-En),2)+4*(y)*pow(Cmn,2)),0.5));
					   
					return Egn_1;
				

}

//**********************************************************************************
//Función para el Cálculo del gap mediante la ley de vegard para un solo valor de x:
//**********************************************************************************

double vegard_law(float EgrBinary_1, float EgrBinary_2, float CTernary_r, double x){
				  
				    double EgrTernary_c;
					EgrTernary_c=EgrBinary_1*x+EgrBinary_2*(1-x)-CTernary_r*x*(1-x);
					
					return EgrTernary_c;
				  
				
}

int main(){

float aSi=5.431; //Lattice parameter for silicon
float aGaAs=5.65325; // lattice parametre a for GaAs
float aGaP=5.4505; // lattice parametre a for GaP
float aGaN=4.52;// lattice parametre a for GaN
double x_prima; // contenido arsenico y fosoforo normalizado
double y_prima;
double y; //contenido en nitrogeno
double x; //contenido en arsenico
double z; //contenido en fosforo	
double VBO; //Band valence offset del GaAsP

double Em_EgrGaAsP; // valor del gap para el GaAsxP1-x con x nomralizado para gamma valley que nos servirá de Em para el BAC
double Eg_GaPAsN;
float En;
float E_nitrogeno=2.18; //ev;
float EgrGaAs=1.424; //valor del Eg para GaAs en gamma valley
float EgrGaP=2.78;//valor del Eg para GaP en gamma valley
float CGaAsPr=0.20;//Bowing parameter en gamma valley

float Cmn_GaAs=2.7;  //todo en eV
float Cmn_GaP=3.05;  //todo en eV
float Cmn;
float VBO_GaAs = -0.8; //Band valence offset del GaAs
float VBO_GaP = -1.27;//Band valence offset del GaP

ofstream f_out("y_x_gap_bulk.dat", ios::out);

for(x=0;x<1.001;x+=0.01){



y=(aSi-aGaP+x*(aGaP-aGaAs))/(aGaN-aGaP);//contenido en nitrogeno lattice match al silicio



x_prima = x/(1-y);	// contenido normalizado de x, denominado x'



cout<<"x: "<<x<<", y: "<<y<<", 1-x-y: "<<1-x-y<<", x': "<<x_prima<<NL;



Em_EgrGaAsP = vegard_law(EgrGaAs, EgrGaP, CGaAsPr, x_prima);// calculo de Em mediante la ley de vegard teniendo como matriz GaAsP con x normalizada




VBO = VBO_GaAs*x_prima + VBO_GaP*(1-x_prima); // calculo de los band offsets del GaAsP ponderizados

En = E_nitrogeno - VBO + VBO_GaP;	//Calculo de En como diferencia de la horizontal En y la recta VBO. Origen de energias en el GaP

cout<< "Em_EgrGaAsP: "<<Em_EgrGaAsP<<", En: "<<En<<NL;


Cmn = Cmn_GaAs*x_prima + Cmn_GaP*(1-x_prima);// Cmn ponderizado


Eg_GaPAsN=BAC_Model(Em_EgrGaAsP, En, Cmn, y); // Calculo del gap del cuaternario a traves del BAC con Em_EgrGaAsP, En y Cmn.

z=1-x-y;	

if(z>0){
	
f_out<<x<<", "<<y<<", "<<Eg_GaPAsN<<NL; // mapa de gaps para un x, y determinado 

}
cout<< "Eg_GaPAsN: "<<Eg_GaPAsN<<NL;

cout<<NL;
}
}

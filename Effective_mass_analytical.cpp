#include <iostream>
#include <vector>
#include <cmath>
using namespace std;
#define NL '\n'
#include <fstream> 

//*****************************************************************************
//Funcion para el Calculo del gap mediante BAC model para un solo valor de y:
//*******************************************************************************
double BAC_Model(float Em, float En, float Cmn, double y){
                 
				    double Egn_1;
					Egn_1=0.5*((Em+En)-pow((pow((Em-En),2)+4*(y)*pow(Cmn,2)),0.5));
					   
					return Egn_1;
				

}

//**********************************************************************************
//Funcion para el Calculo del gap mediante la ley de vegard para un solo valor de x
//**********************************************************************************

double vegard_law(float EgrBinary_1, float EgrBinary_2, float CTernary_r, double x){
				  
				    double EgrTernary_c;
					EgrTernary_c=EgrBinary_1*x+EgrBinary_2*(1-x)-CTernary_r*x*(1-x);
					
					return EgrTernary_c;
				  
				
}

//**********************************************************************************
//Funcion para el Calculo del gap mediante la ley de vegard para un solo valor de x
//**********************************************************************************

double Effective_Mass_N_ter(float m_M, float En, float Em, float Cmn, double y){
				  
				    float eff_mass_N;
				    				    
					eff_mass_N = m_M /(1-(Em-En)/sqrt(pow(Em-En,2)+4*y*pow(Cmn,2)));
					
					return eff_mass_N;
				  
				
}
//**********************************************************************************
//Funcion para el Calculo del gap mediante la ley de vegard para un solo valor de x
//**********************************************************************************

double Effective_Mass_N_qua(double EgrBinary_1, double EgrBinary_2, float CTernary_r, float Cmn_binario_1, float Cmn_binario_2, 
                            float VBO_binario_1, float VBO_binario_2, float M_ternary_r_1, float M_ternary_r_2, float C_M_ternary_r, 
							double x, double y, char *alloy){
				  
				    float eff_mass_N;
				    
					float E_nitrogeno=2.18;
				    
				    double x_prima;
				    
				    double En;
				    double Em;
				    
				    double Cmn;
				    
				    float VBO;
				    
				    double m_M;
				    
				    
					
				    
				       if(alloy =="GaAsPN"){
					    
						   x_prima = x/(1-y);
				       }else{
				    	
				    	   x_prima=x;
					   }
				
					
				    Em = vegard_law(EgrBinary_1, EgrBinary_2,  CTernary_r,  x),
				    
				    
				    VBO = VBO_binario_1*x_prima + VBO_binario_2*(1-x_prima);
				    
				    En = E_nitrogeno - VBO + VBO_binario_2;
				    
				    Cmn = Cmn_binario_1*x_prima + Cmn_binario_2*(1-x_prima);
				    
				    m_M = vegard_law( M_ternary_r_1,  M_ternary_r_2, C_M_ternary_r, x);
				    
				    
					eff_mass_N = m_M /(1-(Em-En)/sqrt(pow(Em-En,2)+4*y*pow(Cmn,2)));
					
					return eff_mass_N;
				  
				
}
int main(){
double x;	
double y;
ofstream f_out("eff_mass.dat", ios::out);// archivo que guarda un la masa effectiva en funcion del contenido del material.
ofstream f_out_1("eff_mass_qua_lattice_match.dat", ios::out);// archivo que guarda un la masa effectiva de un cuaternario
                //lattice match al silicio en funcion de y.
//*********************************************///
//Parametros para GaP1-yNy:
//********************************************///
float Cmn_GaP=3.05;  //todo en eV
float Em_GaPN=2.78; //gap en el centro de la zona de Brillouin. El gap pasa de indirecto (gamma-X) a directo (gamma) cuando se aÃ±ade nitrÃ³geno.
float En_GaPN=2.18;   //referencias
float VBO_GaP = -1.27;//Band valence offset del GaP
float mN_lh_GaP=0.199;// masa de huecos ligeros en la banda de valencia
float mN_hh_GaP=0.326; // masa de huecos pesados en la banda de valencia
float mN_e_GaP=0.13; //masa effectiva del electrón en la banda de conducción

//*********************************************///
//parametros para GaAs1-yNy:
//********************************************///
float Cmn_GaAs=2.7;  //todo en eV
float Em_GaAsN=1.42; //gap en el centro de la zona de Brillouin. El gap pasa de indirecto (gamma-X) a directo (gamma) cuando se aÃ±ade nitrÃ³geno.
float En_GaAsN=1.65;   //referencias
float VBO_GaAs = -0.8; //Band valence offset del GaAs
float mN_lh_GaAs=0.096;// masa de huecos ligeros en la banda de valencia
float mN_hh_GaAs=0.350; // masa de huecos pesados en la banda de valencia
float mN_e_GaAs=0.067; //masa effectiva del electrón en la banda de conducción
//******************************************************///
//parametros para GaAsxP1-x mediante la ley de vegard:
//*****************************************************///
float C_mN_e_GaAsP=0.0086;
float C_mN_lh_GaAsP=0;
float C_mN_hh_GaAsP=0;

// los 3 gaps de los 2 binarios, GaAS y GaP: valles gamma, X y L
float EgrGaAs=1.424; //eV - ioffe
float EgXGaAs=1.90;//eV
float EgLGaAs=1.71;//eV

float EgrGaP=2.78; //ioffe
float EgXGaP=2.26;//eV
float EgLGaP=2.6;//eV
// parametros de curvatura (bowing) de los gaps 
float CGaAsPr=0.20;
float CGaAsPX=0.240;
float CGaAsPL=0.16;

//******************************************************///
//parametros para GaInxP1-x mediante la ley de vegard:
//*****************************************************///
// los 3 gaps de los 2 binarios, GaAS y GaP: gamma X y L

float EgrInP=1.34; //Ioffe: Gap between gamma valley and the heavy holes (or light) band edge for InP in eV
float EgXInP=2.19; //Ioffe: Gap between X valley and the heavy holes (or light) band edge for InP in eV
float EgLInP=1.93; //Ioffe: Gap between L valley and the heavy holes (or light) band edge for InP in eV
float VBO_InP = -0.94;//Band valence offset del InP
float Cmn_InP=3;  //todo en eV
// parametros de curvatura (bowing) de los gaps 

float CGaInPr=0.65;
float CGaInPX=0.2;
float CGaInPL=0.34;	


float C_mN_e_GaInP=0.051;
float C_mN_lh_GaInP=0;
float C_mN_hh_GaInP=0;

float mN_lh_InP=0.121;// masa de huecos ligeros en la banda de valencia
float mN_hh_InP=0.532; // masa de huecos pesados en la banda de valencia
float mN_e_InP=0.0795; //masa effectiva del electrón en la banda de conducción

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


	
   for(x=0 ;x<1.01;x+=0.01){
	
      y=(aSi-aGaP+x*(aGaP-aGaAs))/(aGaN-aGaP) ;
   
      f_out_1<<x<<", "<<y<<", "<<Effective_Mass_N_qua(EgrInP, EgrGaP, CGaInPr, Cmn_GaAs, Cmn_GaP, VBO_InP, VBO_GaAs,
                           mN_e_InP,mN_e_GaP,  C_mN_e_GaInP, x, y, "GaInPN")<<", "<<vegard_law(mN_lh_InP, mN_lh_GaP, C_mN_lh_GaInP, x)<<", "<<
	                         vegard_law(mN_hh_InP, mN_hh_GaP, C_mN_hh_GaInP, x)<<NL;
      //f_out_1<<x<<", "<<y<<", "<<Effective_Mass_N_qua(EgrInP, EgrGaP, CGaInPr, Cmn_GaAs, Cmn_GaP, VBO_InP, VBO_GaAs,
                           //mN_e_InP,mN_e_GaP,  C_mN_e_GaInP, x, y, "GaInPN")<<NL;
                            
      for(y=0 ;y<1.001;y+=0.001){
        
   
   
      cout<<"x: "<<x<<", y: "<<y<<", eff_mass: "<<Effective_Mass_N_qua(EgrInP, EgrGaP, CGaInPr, Cmn_GaAs, Cmn_GaP, VBO_InP, VBO_GaAs,
                           mN_e_InP,mN_e_GaP,  C_mN_e_GaInP, x, y, "GaInPN")<<NL;
    
                                                    
      f_out<<y<<", "<<Effective_Mass_N_ter(mN_e_GaAs, En_GaAsN, Em_GaAsN, Cmn_GaAs, y)<<", "<<mN_lh_GaAs<<", "<<mN_hh_GaAs<<NL;
                            
       
		
	             
   }

}



}

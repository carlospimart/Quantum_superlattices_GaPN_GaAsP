#include <iostream>
#include <vector>
#include <cmath>
#include "vegard_law.cpp"
#include "band_anti-crossing.cpp"
using namespace std;
#define NL '\n'

double Calculator_Band_offset_General(double x, double y, char *type, char *SL_type){


//float x=0.42;//As content in GaAsP
//float y=0.0584466;//%N content in GaPN

float z=0.27;
float Cmn=3.05;  //todo en eV
float Em=2.78; //gap en el centro de la zona de Brillouin. El gap pasa de indirecto (gamma-X) a directo (gamma) cuando se añade nitrógeno.
float En=2.18;   //referencias
float EgrGaAs=1.424; //ioffe
float EgXGaAs=1.90;
float EgLGaAs=1.71;
float EgrGaP=2.78; //ioffe
float EgXGaP=2.26;
float EgLGaP=2.6;
float CGaAsPr=0.20;
float CGaAsPX=0.240;
float CGaAsPL=0.16;

float EgrInP=1.34; //Ioffe: Gap between gamma valley and the heavy holes (or light) band edge for InP in eV
float EgXInP=2.19; //Ioffe: Gap between X valley and the heavy holes (or light) band edge for InP in eV
float EgLInP=1.93; //Ioffe: Gap between L valley and the heavy holes (or light) band edge for InP in eV

float CGaInPr=0.65;
float CGaInPX=0.2;
float CGaInPL=0.34;

float VBO_GaP = -1.27;//Band valence offset del GaP
float VBO_GaAs = -0.8; //Band valence offset del GaAs
float VBO_InP = -0.94;//Band valence offset del InP;

float Cmn_GaP=3.05;  //todo en eV
float Em_GaPN=2.78; //gap en el centro de la zona de Brillouin. El gap pasa de indirecto (gamma-X) a directo (gamma) cuando se aÃ±ade nitrÃ³geno.
float En_GaPN=2.18;   //referencias


float AffGaAs=4.07; // Electronic Affinity for GaAs
float AffGaP=3.8; // Electronic Affinity for GaP
float AffGaN=4.1; // Electronic Affinity for GaN
float mGaAs=0.067;  
float mGaP=0.13;

float Delta_SO_GaP= 0.08;
float Delta_SO_GaAs=0.341;
float Delta_SO_InP=0.108;

float VBO_GaN= -2.64;
float Delta_SO_GaN=-0.017;

float meffGaAsP;
float Egn;
float EgrGaAsP;
float EgXGaAsP;
float AffGaAsP;
float AffGaPN;
float EgMED;
float CBO;
float VBO;
float value;


//meffGaAsP=x*mGaAs+(1-x)*mGaP;
//Egn=0.5*(Em+En-pow(((pow((Em-En),2)+4*(y)*pow(Cmn,2))),0.5));

//EgrGaAsP=EgrGaAs*x+EgrGaP*(1-x)-CGaAsPr*x*(1-x);
//EgXGaAsP=EgXGaAs*x+EgXGaP*(1-x)-CGaAsPX*x*(1-x);

//AffGaAsP=AffGaAs*x+AffGaP*(1-x); //Electronic Affinity for GaAsP
//AffGaPN=AffGaP*(1-y)+AffGaN*(y); //Electronic Affinity for GaPN

//EgMED=z*EgXGaAsP+(1-z)*Egn;

//CBO=AffGaAsP - AffGaPN;// Difference between Electronic Affinities that makes Conduction band offset for GaPN/GaAsP

//VBO=(AffGaAsP+EgrGaAsP) - (AffGaPN+Egn);// Difference between Electronic Affinities that makes Conduction band offset for GaAsP/GaPN
if(SL_type == "GaAsP/GaPN"){

    
	VBO = VBO_GaP*(1-x) + VBO_GaAs*x  - (VBO_GaP*(1-y) +VBO_GaN*y);//1/3.0*(Delta_SO_GaP*(1-x)+Delta_SO_GaAs*x) - 
    
	//(VBO_GaP*(1-y) + VBO_GaN*y + 1/3.0*(Delta_SO_GaP*(1-y)+Delta_SO_GaN*y)) ;//Valence Band offset

    CBO= vegard_law_func(EgrGaAs, EgrGaP, CGaAsPr, x) - BAC_Model_func(Em, En, Cmn_GaP, y) + VBO ;//conduction band offset
    
     //cout<<"VBO: "<<VBO<<", "<<"CBO: "<<CBO<<", GaAsP: "<<vegard_law_func(EgrGaAs, EgrGaP, CGaAsPr, x)<<"GaPN: "<< BAC_Model_func(Em, En, Cmn_GaP, y)<<NL;

}else if(SL_type == "GaInP/GaPN"){
	
    VBO = VBO_GaP*(1-x) + VBO_InP*x - (VBO_GaP*(1-y) + VBO_GaN*y); //+ 1/3.0*(Delta_SO_GaP*(1-x)+Delta_SO_InP*x) - 
	 //(VBO_GaP*(1-y) + VBO_GaN*y + 1/3.0*(Delta_SO_GaP*(1-y)+Delta_SO_GaN*y)) ;//Valence Band offset

    CBO=vegard_law_func(EgrInP, EgrGaP, CGaInPr, x) - BAC_Model_func(Em, En, Cmn_GaP, y) + VBO ; //conduction band offset
}

if(type == "CBO"){
	value=CBO;
}else if(type == "VBO"){
    value=VBO; 
}else{
	value=NAN; 
}

return value;

}


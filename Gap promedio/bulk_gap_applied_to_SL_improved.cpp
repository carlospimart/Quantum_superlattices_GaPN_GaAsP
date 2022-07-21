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
//Funcion para el Calculo del gap mediante la ley de vegard para un solo valor de x:#
//**********************************************************************************

double vegard_law(float EgrBinary_1, float EgrBinary_2, float CTernary_r, double x){
				  
				    double EgrTernary_c;
					EgrTernary_c=EgrBinary_1*x+EgrBinary_2*(1-x)-CTernary_r*x*(1-x);
					
					return EgrTernary_c;
				  
				
}

//*********************************************************************************
//Funcion para el Calculo del espesor relativo de una superred en funcion de x e y 
//*********************************************************************************
double REL_Func(float a_binario1_a, float a_binario1_b, float a_binario2_a, float a_binario2_b, double x, double y){
                 
				    double REL;
				    float aSi=5.431;
				    REL=((aSi-(a_binario2_a*(1-y)+y*a_binario2_b)))/((a_binario1_a*x+(1-x)*a_binario1_b)-(a_binario2_a*(1-y)+y*a_binario2_b));
					   
					return REL;
				

}
//****************************************************************************************************
//Funcion para el Calculo del gap de una superred mediante el bulk pesado por sus espesores relativos   
//****************************************************************************************************
double SL_Gap_bulk(float a_binario1_a, float a_binario1_b, float a_binario2_a, float a_binario2_b, double x, double y,
                   float Egr_binario1, float Egr_binario2, float C_Ternario_r, float VBO_binario_1, float VBO_binario_2, 
				   float Cmn_binario_1, float Cmn_binario_2, float rel_thick, bool n, char *alloy){
                 
                    double gap; // variable donde guardamos el valor del gap
                    float Cmn; //elemento matriz
                    float VBO; //variable para guardar el band offset 
                    float En; //Nivel loalizado del N
                    float REL; // variable para guardar el espesor relativo a través de la función REL_Func 
                    float x_rel; //variable de composición x pesada por REL
                    float y_rel; //variable de composición y pesada por REL
                    float x_prima; // variable x 
                    float E_nitrogeno=2.18; //nivel del N para x=0
                    
                    double Em_Egr_ternario; // variable donde calculamos el valor Em del ternario matriz
                    
				    if(n==true){//condición para decidir si queremos calcular el espesor relativo para un x e y dados, o elegimos 
					            //uno impuesto por nosotros
					    REL=REL_Func(a_binario1_a, a_binario1_b, a_binario2_a, a_binario2_b, x, y);
                    }else{
                    	REL =rel_thick;
					}	
					
                    x_rel=x*REL; // contenido en Arsenico pesado con el espesor relativo de GaAsP
                    y_rel=y*(1-REL);// contenido en Nitrogeno pesado con el espesor relativo de GaPN
		
					
					if(alloy =="GaAsPN"){ // condición para decidir si normalizamos o no; si se tratara de GaInPN NO HABRÍA NECESIDAD DE NORMALIZAR
					    
						x_prima = x_rel/(1-y_rel);
				    }else{
				    	
				    	x_prima=x_rel;
					}
							
				    Em_Egr_ternario = vegard_law(Egr_binario1, Egr_binario2, C_Ternario_r, x_prima);
				    
				    VBO = VBO_binario_1*x_prima + VBO_binario_2*(1-x_prima);
				    
				    En = E_nitrogeno - VBO + VBO_binario_2;
				    cout<<"En: "<<En<<NL;
				    Cmn = Cmn_binario_1*x_prima + Cmn_binario_2*(1-x_prima);
				    
				    gap=BAC_Model(Em_Egr_ternario, En, Cmn, y_rel);
				    
					
					return gap;
				

}
//*******************************************************************************************************************************
//Funcion para el Calculo del gap de una superred mediante el bulk pesado por sus espesores relativos para el caso de GaAsN/GaPN   
//*******************************************************************************************************************************
double SL_Gap_bulk_2(float a_binario1_a, float a_binario1_b, float a_binario2_a, float a_binario2_b, double x, double y,
                   float Em_1, float En_1, float Em_2, float En_2, float Cmn_1, float Cmn_2, float C_Ternario_r, 
				   float rel_thick, bool n){
				   
				float x_rel;
                float y_rel;
                float x_prima;
                
                double Eg_1_N;
                double Eg_2_N;
                double REL;
                double gap;
                
				if(n==true){//condiciÃ³n para decidir si queremos calcular el espesor relativo para un x e y dados, o elegimos 
					            //uno impuesto por nosotros
					    REL=REL_Func(a_binario1_a, a_binario1_b, a_binario2_a, a_binario2_b, x, y);
                    }else{
                    	REL =rel_thick;
					}	
					
                    x_rel=x*REL; // contenido en Arsenico pesado con el espesor relativo de GaAsP
                    y_rel=y*(1-REL);// contenido en Nitrogeno pesado con el espesor relativo de GaPN
                    
				Eg_1_N=BAC_Model(Em_1, En_1, Cmn_1, x_rel);
				    
				Eg_2_N=BAC_Model(Em_2, En_2, Cmn_2, y_rel);
					
				x_prima = x_rel + x_rel*y/(1-x_rel); //NORMALIZACIÓN DE LA variable  x
					
				gap = vegard_law(Eg_1_N, Eg_2_N, C_Ternario_r, x_prima);
}
int main(){
//*********************************************///
//Parametros para GaP1-yNy:
//********************************************///
float Cmn_GaP=3.05;  //todo en eV
float VBO_GaP = -1.27;//Band valence offset del GaP
float Em_GaPN=2.78; //gap en el centro de la zona de Brillouin. El gap pasa de indirecto (gamma-X) a directo (gamma) cuando se aÃ±ade nitrÃ³geno.
float En_GaPN=2.18;   //referencias
//double Egn_GaPN[dimy_1];//EnergÃ­a del gap de la banda BAC E- del GaPN del nitrogeno
//double Egm_GaPN[dimy_1];//EnergÃ­a del gap de la banda BAC E- del GaPN del nitrogeno


//*********************************************///
//parametros para GaAs1-yNy:
//********************************************///
float Cmn_GaAs=2.7;  //todo en eV
float Em_GaAsN=1.42; //gap en el centro de la zona de Brillouin. El gap pasa de indirecto (gamma-X) a directo (gamma) cuando se aÃ±ade nitrÃ³geno.
float En_GaAsN=1.65;   //referencias
float VBO_GaAs = -0.8; //Band valence offset del GaAs
//double Egn_GaAsN[dimx_1]; //EnergÃ­a del gap de la banda BAC E- del GaAsN del nitrogeno
//double Egm_GaAsN[dimx_1];// EnergÃ­a del gap de la banda BAC E+ del GaAsN del nitrogeno
//******************************************************///
// Parametros para GaAsxP1-x mediante la ley de vegard:
//*****************************************************///

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
//Calculo del gap de GaInxP1-x mediante la ley de vegard:
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


//*******************************************************************//
//*********************superlattice gap promedio******************************//
//*******************************************************************//



ofstream f_out("y_x_gap_SL.dat", ios::out);// archivo que guarda un mapa de gaps para un valor de x e y dados
ofstream f_out_1("rel_thick_for_gap_SL.dat", ios::out);//archivo declarado para guardar una recta de x e y para un espesor relativo constante
ofstream f_out_2("x_y_thickness_SL.dat", ios::out);//archivo declarado para guardar un mapa de espesores relativos para un valor de x e y dados
ofstream f_out_4("y_x_SL_gap.dat", ios::out);// archivo declarado para grabar x,y e el espesor realtivo para un gap dado
ofstream f_out_5("x_ter_gap.dat", ios::out);// archivo declarado para grabar x o y y el gap correspondiente de un ternario

float SL_Gap0;// calculo del gap para x0, y0
float rel_thick=0.25;// valor del espesor relativo para un x e y dados
float SL_Gapf;// valor del gap para xf, yf
float SL_Gapm;//velor intermedio entre sp_gap1 y sa_gap2
double x;
double y;
double x0=0;
double xf=0;
double yf=0;
double y0=0;
double xm=0;
double ym=0;
double y_for_rel;
double Eg_GaPAsN;
double Eg_GaInPN;
double Eg_GaPN;
double temp = 0;


float SL_Gap_promedio;
float SL_GAP_diferencia;

float aSi=5.431; //Amstrongs - Lattice parameter for silicon
float aGaAs=5.65325; //Amstrongs - lattice parametre a for GaAs
float aGaP=5.4505; //Amstrongs - lattice parametre a for GaP
float aGaN=4.52;//Amstrongs -  lattice parametre a for GaN
float aInP=5.8687;//Amstrongs -  lattice parametre a for InP. Ioffe
float aInN=5.004;//Amstrongs -  lattice parametre a for InN
float aGaIn=5.004;//Amstrongs -  lattice parametre a for GaIn

float SL_Gap=1.9;

cout<<"Do you want to calculate the gaps for a nitride ternary in a range of composition?"<<NL;
cout<<"Type 1 if you want to carry on with this process or if you don't want to."<<NL;
cout<<"Type 2 if you don't want to."<<NL;

int input_ter;

cin>>input_ter;

switch(input_ter){// switch para decidir si  queremos calcular los gaps de algun ternario 

case 1:
    
  for(y=0;y<1;y+=0.001){
	
	f_out_5<<y<<", "<<BAC_Model(Em_GaPN, En_GaPN, Cmn_GaP, y)<<NL;
	
  }
break;

case 2:
	cout<<"Let's continuous"<<NL;
	break;
default:
	cout<<"wrong entrance"<<NL;
}

int input_1;
cout<<NL;
cout<<"Do you want to create a map of gaps for each x e y or just calculate the gap only for one value of x,y and REL?"<<NL;
cout<<"- Type 1 if you want just to calculate the gap only for one value of x, y and relative thickness."<<NL;
cout<<"- Type 2 if you want to create a map of gaps for each x e y."<<NL;
cout<<"- Type 3 if don't want to carry on with any of these process'"<<NL;

cin>>input_1;

switch(input_1){// switch para decidir si  queremos calcular el mapa de gaps para cada valor x, y 

case 1:
	 cout<<"Introduce the value of x"<<NL;
	 cin>>x;
	
	 cout<<"Introduce the value of relative thickness"<<NL;
	 
	 cin>>rel_thick;
	 
	 y_for_rel=((aSi-rel_thick*(aGaAs*x+(1-x)*aGaP))/((1-rel_thick)*(aGaN-aGaP))) - (aGaP/(aGaN-aGaP)); // recta x, y para un REL constante para GaAsPN.
	 
	 Eg_GaPAsN =  SL_Gap_bulk(aGaAs, aGaP, aGaP, aGaN, x, y_for_rel, EgrGaAs, EgrGaP, CGaAsPr, VBO_GaAs, VBO_GaP, Cmn_GaAs, Cmn_GaP, rel_thick, false,"GaAsPN");
     
	 Eg_GaPN = BAC_Model(Em_GaPN, En_GaPN, Cmn_GaP, y_for_rel);
	 
	 cout<<"x: "<<x<<", "<<", y: "<<y_for_rel<<", "<<", REl: "<<rel_thick<< ", Eg_GaPN: "<<Eg_GaPN <<", Eg_GaPAsN: "<<Eg_GaPAsN<<NL;
	 
break;

case 2:
    for(x=0;x<0.61;x+=0.01){
	
//y_for_rel=((aSi-rel_thick*(aGaAs*x+(1-x)*aGaP))/((1-rel_thick)*(aGaN-aGaP))) - (aGaP/(aGaN-aGaP));//Calculo de las rectas de REL constantes para x e y para el GaInP/GaPN
        y_for_rel=((aSi-rel_thick*(aInP*x+(1-x)*aGaP))/((1-rel_thick)*(aGaN-aGaP))) - (aGaP/(aGaN-aGaP)); //Calculo de las rectas de REL constantes para x e y
        
		Eg_GaInPN =  SL_Gap_bulk(aInP, aGaP, aGaP, aGaN, x, y_for_rel, EgrInP, EgrGaP, CGaInPr, VBO_InP, VBO_GaP, Cmn_InP, Cmn_GaP, 0.5, false,"GaInPN");
        
		f_out_1<<x<<", "<<y_for_rel<<", "<<Eg_GaInPN<<NL;
   
   for(y=0.02 ;y<0.151;y+=0.001){
 	
      f_out_2<<x<<", "<<y<<", "<<REL_Func(aGaAs, aGaP, aGaP, aGaN, x, y)<<NL; 
	
      //Eg_GaPAsN =  SL_Gap_bulk(aGaAs, aGaP, aGaP, aGaN, x, y, EgrGaAs, EgrGaP, CGaAsPr, VBO_GaAs, VBO_GaP, Cmn_GaAs, Cmn_GaP, NAN, true,"GaAsPN");
    
    Eg_GaInPN =  SL_Gap_bulk(aInP, aGaP, aGaP, aGaN, x, y, EgrInP, EgrGaP, CGaInPr, VBO_InP, VBO_GaP, Cmn_InP, Cmn_GaP, NAN, true,"GaInPN"); 
	   
	//Eg_GaPAsN =  SL_Gap_bulk_2(aGaAs, aGaN, aGaP, aGaN, x, y, Em_GaAsN, En_GaAsN, Em_GaPN, En_GaPN, Cmn_GaAs, Cmn_GaP, CGaAsPr, NAN, true); 
	
		       
       f_out<<x<<", "<<y<<", "<<Eg_GaInPN<<NL;

       cout<< "Eg_GaInPN: "<< Eg_GaInPN<<NL;

	cout<<NL;
}
}
break;
case 3:
		cout<<"the process was finished"<<NL;
	break;
default: 
		    cout<<"Please, type a selection correctly"<<NL;
	
}
int input;
cout<<NL;
cout<<"Do you want to calculate the values of x, y and relative thcikness for a given value of band-gap energy?"<<NL;
cout<<"- Type 1 if want to carry on with this process."<<NL;
cout<<"- Type 2 if don't want to'"<<NL;
cin>>input;

switch(input){// switch para decidir si ademas queremos calcular con el metodo numerico la cuva de gap constante para cada valor x, y y Rel 

case 1: 
     
	 for(rel_thick=0; rel_thick<1;rel_thick+=0.01){
	
     cout<<"rel_thick: "<<rel_thick<<NL;	
     x0=0;// contenido en Arsenico del GaAsP
     xf=0.6;
     temp=1000;
     // realthick es el porcentaje de la capa de arsenico GaAsP en la superred 
     y0=((aSi-rel_thick*(aGaAs*x0+(1-x0)*aGaP))/((1-rel_thick)*(aGaN-aGaP))) - (aGaP/(aGaN-aGaP));

     yf=((aSi-rel_thick*(aGaAs*xf+(1-xf)*aGaP))/((1-rel_thick)*(aGaN-aGaP))) - (aGaP/(aGaN-aGaP));


     SL_Gap0= SL_Gap_bulk(NAN, NAN, NAN, NAN, x0, y0, EgrGaAs, EgrGaP, CGaAsPr, VBO_GaAs, VBO_GaP, Cmn_GaAs, Cmn_GaP, rel_thick, false,"GaAsPN");          
     cout<<"SL_Gap0: "<< SL_Gap0<<NL;
     SL_Gapf= SL_Gap_bulk(NAN, NAN, NAN, NAN, xf, yf, EgrGaAs, EgrGaP, CGaAsPr, VBO_GaAs, VBO_GaP, Cmn_GaAs, Cmn_GaP, rel_thick, false,"GaAsPN"); 
     cout<<"SL_Gapf: "<< SL_Gapf<<NL;

     SL_Gap_promedio=(SL_Gap0+SL_Gapf)/2;

     SL_GAP_diferencia=abs((SL_Gap0-SL_Gapf)/2);


     if(SL_Gap_promedio-SL_GAP_diferencia<SL_Gap   &&   SL_Gap<SL_Gap_promedio+SL_GAP_diferencia){
   
    
     while (abs(temp)>10E-6){
	
     y0=((aSi-rel_thick*(aGaAs*x0+(1-x0)*aGaP))/((1-rel_thick)*(aGaN-aGaP))) - (aGaP/(aGaN-aGaP));


     yf=((aSi-rel_thick*(aGaAs*xf+(1-xf)*aGaP))/((1-rel_thick)*(aGaN-aGaP))) - (aGaP/(aGaN-aGaP));

     xm=(x0+xf)/2;

     ym=((aSi-rel_thick*(aGaAs*xm+(1-xm)*aGaP))/((1-rel_thick)*(aGaN-aGaP))) - (aGaP/(aGaN-aGaP));

     //cout<<"GaPN0: "<<BAC_Model(Em_GaPN, En_GaPN, Cmn_GaP, y0)<<", GaAsP0: "<<vegard_law(EgrGaAs, EgrGaP, CGaAsPr, x0)<<NL;
 
     SL_Gap0= SL_Gap_bulk(NAN, NAN, NAN, NAN, x0, y0, EgrGaAs, EgrGaP, CGaAsPr, VBO_GaAs, VBO_GaP, Cmn_GaAs, Cmn_GaP, rel_thick, false,"GaAsPN");
         
     //cout<<"GaPNf: "<<BAC_Model(Em_GaPN, En_GaPN, Cmn_GaP, yf)<<", GaAsPf: "<<vegard_law(EgrGaAs, EgrGaP, CGaAsPr, xf)<<NL;     
   
     SL_Gapf= SL_Gap_bulk(NAN, NAN, NAN, NAN, xf, yf, EgrGaAs, EgrGaP, CGaAsPr, VBO_GaAs, VBO_GaP, Cmn_GaAs, Cmn_GaP, rel_thick, false,"GaAsPN");
         
     SL_Gapm= SL_Gap_bulk(NAN, NAN, NAN, NAN, xm, ym, EgrGaAs, EgrGaP, CGaAsPr, VBO_GaAs, VBO_GaP, Cmn_GaAs, Cmn_GaP, rel_thick, false,"GaAsPN");
 
     if(SL_Gap0>SL_Gapf){
	
	 if(SL_Gap>SL_Gapm){
		
		  xf=xm;
		
	   
	  }else{
	  	 x0=xm;
	  }
	
     }else{
	
	 if(SL_Gap>SL_Gapm){
		
		  x0=xm;
		
	   
	 }else{
	   	  xf=xm;
	 }
    }
       temp=SL_Gap0-SL_Gapf;
   
    cout<<"temp: "<<temp<<NL;
   
    }

    if(((y0+yf)/2)<=1){// condición impuesta para descartar valores de y mayores que 1

    cout<<"xdef: "<<(x0+xf)/2<<", ydef: "<<(y0+yf)/2<<", SL_Gap0: "<<SL_Gap0<<NL;
    f_out_4<<(x0+xf)/2<<", "<<(y0+yf)/2<<", "<<rel_thick<<NL;
    }else{
    cout<<"No solutions"<<NL;	
    }

    }else{
	
    cout<<"out of range"<<NL;
	
    }
    }
break;

case 2:
	
	cout<<"the process was finished"<<NL;
	break;
default: 
		    cout<<"Please, type a selection correctly"<<NL;
}

}

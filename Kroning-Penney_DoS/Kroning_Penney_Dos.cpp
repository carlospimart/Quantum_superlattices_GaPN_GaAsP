#include <iostream>
#include <cmath>
#include <math.h>
#include <fstream>
#include "Calculator_Band_offset_General.cpp"
#include "protogaps_function.cpp"
#include "REL_func.cpp"
//#include "Band_Offset_GaAsN_GaPN.cpp"
#include "Effective_mass_analytical_function.cpp"
using namespace std;
#define NL '\n'


int main(){

int hole_type=3;// parameters to decide if we want to compute the heavy hole energy band (hole_type = 3) or light hole energy band (hole_type = 2)
//*************************************************************//
//Effective mass: 
//
//Effective_mass_analytical_function(double x, char *type, int input)
//parameter "input" = 1: Electron effective mass 
//parameter "input" = 2: light hole effective mass
//parameter "input" = 3: heavy hole effective mass 
//*************************************************************//
char *type = "GaAsP/GaPN"; //parameter to decide what is the superlattice that we want to simulate. 
// The superlattice availiable in the system is:
// - GaAsP/GaPN
// - GaInP/GaPN
double x=0.17; //Com 
double y=0.04; //composición en nitrogeno
double REL=0.25;// Espesor relativo de la superred
double d=4.9E-9; //anchura de un periodo (pozo+barrera) (recomendable hasta 17 nanometros)
double a=d*(1-REL);  //anchura del pozo en nanometros

double vector;
char *type_2; // variable para indicar en forma de parametro el smiconductor en cuestión. se utiliza para el calculo de la masa efectiva
              // en la función "Effective_mass_analytical_function.cpp"
double delta_gap; //variable utilizada para acumular el valor E0 en k=0 en  la banda de conducción y de valencia
double gap;// variable del gap de la heteroestructura (delta_gap + protogap)
int NumBands=3; //numero de minibandas a calcular
double Vo; //variable para  guardar los Band offset, estos son obtenidos a traves de la función: "Calculator_Band_offset_General.cpp"
float mb;// masa efectiva en la barrera (eV)
float mw;// masa efectiva en el pozo (eV)
double me=9.11E-31; //masa efectiva del electón en kg 

int l=0;
int rest=0;

double hbar=1.05457E-34;   //J.s
// Existen soluciones a la ecuacion de onda cuando la siguiente ecuacion
// tiene solucion (relación de dispersión)
int i=0;
int t;
double E=0; // variable para acumular los niveles de energías delas minibandas
int n;
double k;
const double pi = acos(-1); //valor de Pi
double alpha; //Kx para V0 <0
double kappa;//KAPPAx para V0 < 0
double beta;//KAPPAx para V0>0
ofstream f_out_1("KP_k1_E1__CB.dat", ios::out);
ofstream f_out_2("KP_k1_E1__VB.dat", ios::out);
ofstream f_out_3("KP_DOS_CB.dat", ios::out);
ofstream f_out_4("KP_DOS_VB.dat", ios::out);
ofstream f_out_5("MAP_GAP_KP.dat", ios::out);
double delta_k = (pi/d)/100; // numero de pasos arbitrario de 0 a pi/d (primera zona de broulline). En este caso 100 pasos



int mdim=NumBands*100 + 1; //dimensión de las marices K, E1, E_VB Y E_CB

int m;
cout<<"mdim: "<<mdim<<NL;
//double** K = new double*[mdim];
//double** E1 = new double*[mdim];

cout<<"delta_k: "<<delta_k<<NL;
double K[mdim]; // MATRIZ QUE GUARDA EL VALOR DEL vector de la red reciproca
double E1[mdim];// MATRIZ QUE GUARDA EL VALOR DE las energia 
double E_VB[mdim];// MATRIZ QUE GUARDA EL VALOR DEL las energías para VB
double E_CB[mdim];// MATRIZ QUE GUARDA EL VALOR DEL las energías para CB
double DOS[2*mdim];// MATRIZ QUE GUARDA EL VALOR DEL la densidad de estados
double deriv[mdim];// MATRIZ QUE GUARDA EL VALOR DEL la derivada de las energías con respecto al momento cristalino



if(type=="GaAsP/GaPN"){//codigo para  decidir si type_2 es igual a "GaAsP" o "GaInP" segun el tipo de superred que tengamos
	
	type_2 = "GaAsP";
	
}else if(type=="GaInP/GaPN"){

    type_2 = "GaInP";
}

for(m=1;m<3;m+=1){
cout<<"m: "<<m<<NL;
i=0;
E=0;
vector = 0;   
    
    

if(m==1){// condicional para diferenciar el calculo de CBO (m=1) o VBO (m=2)
	
	Vo= Calculator_Band_offset_General(x, y, "CBO", type)*1.6E-19; //función que calcula la altura de la barrera (band offset) para el CBO
	
	cout<<"CBO: "<<Calculator_Band_offset_General(x, y, "CBO", type)<<NL;
	
	
	
    
    
    
    
    mb = Effective_mass_analytical_function(x, type_2, 1);//funcion que calculan LA MASA EFECTIVA de la barrera
    
    mw = Effective_mass_analytical_function(y, "GaPN", 1); //funcion que calculan LA MASA EFECTIVA del pozo 
    
}else{
	Vo=Calculator_Band_offset_General(x, y, "VBO", type)*1.6E-19;
	
	
	cout<<"VBO: "<<Calculator_Band_offset_General(x, y, "VBO", type)<<NL; //función que calcula la altura de la barrera (band offset) para el VBO
	
    mw = Effective_mass_analytical_function(x, type_2, hole_type);//funcion que calculan LA MASA EFECTIVA del pozo
    
    mb= Effective_mass_analytical_function(y, "GaPN", hole_type);//funcion que calculan LA MASA EFECTIVA de la barrera
   
    a= d*(REL);
    
    
    
}
cout<<"mw: "<<mw<<NL;
    
cout<<"mb: "<<mb<<NL;
    
cout<<"a: "<<a<<NL;
mw=mw*me; //mw en Kg
mb=mb*me; // mb en Kg

for(n=1;n<NumBands+1;n++){
	cout<<"n: "<<n<<NL;
	cout<<NL;
	for(k=((n-1)*(pi/d));k<(n*(pi/d));k+=delta_k){
		
		//k=0;*/
		
		double temp=1000;
		if(E<=Vo){ // relación de dispersión de eenergías obtenida con el modelo de Kronig Penney
			while (abs(temp)>10E-5){
			E=E+(0.00000001*1.6E-19);	
            alpha=sqrt(2*mw*E/pow(hbar,2));
			kappa=sqrt(2*mb*(Vo-E)/pow(hbar,2));
			temp=cosh(kappa*(d-a))*cos(alpha*a)+((pow(kappa,2)-pow(alpha,2))*sinh(kappa*(d-a))*sin(alpha*(a))/(2*alpha*kappa))-cos(k*d);
			//cout<<"temp: "<<temp<<NL;
			 
			}
			
		}
		else{
			
			while (abs(temp)>10E-5){
            alpha=sqrt(2*mw*E/pow(hbar,2));
            beta=sqrt(2*mb*(E-Vo)/pow(hbar,2));
			temp=cos(beta*(d-a))*cos(alpha*a)-((pow(alpha,2)+pow(beta,2))*sin(beta*(d-a))*sin(alpha*a)/(2*alpha*beta))-cos(k*d);
		    //cout<<"temp: "<<temp<<NL;
			E=E+(0.00000001*1.6E-19);
			
			}
		
            
		}
        
		
		
		
		K[i]=k*d/pi; 
		
		
       E1[i]=E/(1.6E-19); //E en eV
        
        		
	     
	/* codigo para guardar la dispersión de energía del la VB Y CB  en dos archivos diferentes.
	   por alguna razón desconocida afecta al calculo del gap. No tilizar sin revisar y slucionar tal problema
	    if(m==1){//for ploting
			
	       vector +=pow(-1,n-1) * delta_k*(d/pi);
	       
	       E_CB[i]=E/(1.6E-19);
	       
       	   f_out_1<<vector-0.01<<", "<<E/(1.6E-19)<<NL;
       	   
	    }else{
	    	
	    	vector +=pow(-1,n-1) * delta_k*(d/pi);
	    	
	       E_VB[i]=E/(1.6E-19);
	       
       	    f_out_2<<vector-0.01<<", "<<-E/(1.6E-19)<<NL;
		}
	*/
		
		cout<<"i: "<<i<<", k: "<<K[i]<<", E: "<<E1[i]<<", E0: "<<E1[0]<<NL;
        
        i+=1;    
	}
         
cout<<"E0: "<<E1[0]<<NL;    
	}// Numbands
//*************************************************/
//**************DENSITY OF STATES *****************//
//*************************************************//

double increment=0.01;//delta_k*(pi/d); //incremento de k usado en la rel. de dispersión. Se usa para hallar la deriv. numérica

//cout<<"increment: "<<increment<<NL;

//Calculo de la derivada numerica y de la DoS

  deriv[0]=((E1[1]-E1[0]))/increment;
  
    
  for (int i=1;i<mdim-1;i++){
    	
   deriv[i]=((E1[i+1]-E1[i-1])/2)/increment;
    
	}
	
	deriv[mdim-1]=(E1[mdim-1]-E1[mdim-2])/increment;
	
for (int i=0;i<mdim;i++){
    	
   //cout<<i<<" "<<deriv[i]<<NL;
    
}
for(i=0;i<mdim;i++){
	
	

    DOS[l]=(2/pi)*(1/deriv[i]);
    
	l+=1;
    
   if(DOS[l]==INFINITY){
	
       DOS[l]=0;
    }
    
    
   
   
}  	


cout<<"E0: "<<E1[0]<<NL;

delta_gap+=E1[0]; //  acumulación de E0 en delta_gap 
	
cout<<"E(0): "<<E1[0]<<NL;

}// m=values for m 

gap = delta_gap +  protogaps_function(x, y, type); // gap: suma de los E0 mas el protogap 
                                                   //(el cual es obtenido a traves de la llamada a la función protogaps_function(x, y, type))

cout<<"delta_gap: "<<delta_gap<<NL;
cout<<"gap: "<<gap<<NL;

cout<<"x: "<<x<<", y: "<<y<<NL;

//f_out_5<<x<<", "<<y<<", "<<gap<<NL;

//}
//}

for(i=0;i<2*mdim;i++){// codigo donde guardamos el valor de la DoS, E para CB y VB. 
  cout<<"i: "<<i<<NL;  
	if(i>=mdim){
    	rest=mdim;
    	f_out_4<<DOS[i]<<", "<<-E_VB[i-rest]<<NL;
		cout<<"DoS: "<<DOS[i]<<", EV: "<<-E_VB[i-rest]<<NL;
	}else{
        f_out_3<<DOS[i]<<", "<<E_CB[i-rest]+gap<<NL;
        cout<<"DoS: "<<DOS[i]<<", EC: "<<E_CB[i-rest]+gap<<NL;
	}
}




}




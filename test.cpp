#include <iostream>
#include <vector>
#include <cmath>
//#include "band_anti-crossing.cpp"
#include "Calculator_Band_offset_General.cpp"
//#include "protogaps_function.cpp"
using namespace std;
#define NL '\n'
#include <fstream> 

int main(){
	
	cout<<Calculator_Band_offset_General(0.31, 0.051, "VBO","GaAsP/GaPN")<<NL;
	
	//cout<<protogaps_function(0.4, 0.047)<<NL;
	//cout<<Band_Offset_GaInP_GaPN(1, 1, "VBO")<<NL;
	
	//cout<<vegard_law_func(1.34, 2.78, 0.65, 1)<<NL;
	
	//cout<<	BAC_Model_func(2.78, 2.18,3.05, 0.653)<<NL;
	
	//cout<<	BAC_Model_func(1.424, 1.65, 2.7, 0)<<NL;
}

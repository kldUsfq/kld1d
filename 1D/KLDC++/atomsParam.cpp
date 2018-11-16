/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   atomsParam.cpp
 * Author: vrodriguez
 * 
 * Created on September 17, 2018, 10:01 PM
 */

#include "atomsParam.h"
#include <string>

const std::string symbols[32] ={"H","Li","Be","B","C","N","O","F","Na","Mg","Al","Si",
			  "P","S","Cl","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co",
			  "Ni","Cu","Zn","Ga","Ge","As","Se","Br"};
const double bsr [32] = {0.661404096e0,1.370051342e0,0.992106144e0,0.803133545e0,
		   0.661404096e0,0.614160946e0,0.566917797e0,0.472431497e0,  
                   1.700753390e0,1.417294491e0,1.181078743e0,1.039349294e0,  
                   0.944862994e0,0.944862994e0,0.944862994e0,2.078698587e0,  
                   1.700753390e0,3.023568000e0,2.645633000e0,2.551135500e0,  
                   2.645633000e0,2.645633000e0,2.645633000e0,2.551135500e0,  
                   2.551135500e0,2.551135500e0,2.551135500e0,1.228321893e0,  
                   1.138559908e0,1.086592443e0,1.086592443e0,1.086592443e0};
 const int n_bc [32] = {20,25,25,25,25,25,25,25,30,30,30,30,30,30,30,35,35,35,35, 
                  35,35,35,35,35,35,35,35,35,35,35,35,35};


 /**
  * returns the symbol according to the number
  * @param symbolNumber
  * @return 
  */
 std::string getSymbol (int symbolNumber){
	if (symbolNumber>=0 && symbolNumber<32){
		return symbols[symbolNumber];
	}else {
		return " ";
	} 
 }
 
 /**
 * returns the nbc for symbolNumber
 * @param asciiCode of the symbol
 * @return 
 */
int getNbc (int symbolNumber){
	if (symbolNumber>=0 && symbolNumber<32){
		return n_bc[symbolNumber];
	}else {
		return 0;
	}
};
/**
 * returns the  bsr of symbolNumber
 * @param symbolNumber
 * @return 
 */
double getBsr (int symbolNumber){
	if (symbolNumber>=0 && symbolNumber<32){
		return bsr[symbolNumber];
	}else {
		return 0.0;
	}
}
/**
 * finds the symbol and returns the position
 * @param symbol
 * @return position
 */
int getSymbolNumber (std::string symbol){
	int symbolNumber = -1;
	for (int i=0;i<32;i++){
		if (symbol == symbols[i]){
			symbolNumber = i;
			break;
		}
	}
	return symbolNumber;
} 

 


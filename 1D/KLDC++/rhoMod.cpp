/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   rhoMod.cpp
 * Author: vrodriguez
 * 
 * Created on September 24, 2018, 3:57 PM
 */

#include <iostream>
#include "rhoMod.h"
#include <math.h>
#include <string>
#include <vector>


/**
 * gaussian is a calculation of the gaussian method
 * @param NX
 * @param xx
 * @param zz
 * @param ALPHA
 * @return 
 */
long double gaussian (int NX, long double xx, long double yy,long double zz, double ALPHA){
	long double result = 0.0;
	long double product=0.0;
	long double ri=0.0;
	
	
	switch (NX){
	case 1:
		product = 1.0;
		break;
	case 2:
		product = xx;
		break;
	case 3:
		product = yy;
		break;
	case 4:
		product = zz;
		break;
	case 5:
		product = powl(xx,2);
		break;
	case 6:
		product = powl(yy,2);
		break;
	case 7:
		product = powl(zz,2);
		break;
	case 8:
		product = xx*yy;
		break;
	case 9:
		product = xx*zz;
		break;
	case 10:
		product = yy*zz;
		break;
	case 11:
		product = pow(xx,3);
		break;
	case 12:
		product = powl(yy,3);
		break;
	case 13:
		product = pow(zz,3);
		break;
	case 14:
		product = powl(xx,2)*yy;
		break;
	case 15:
		product = powl(xx,2)*zz;
		break;
	case 16:
		product = powl(yy,2)*zz;
		break;
	case 17:
		product = powl(yy,2)*xx;
		break;
	case 18:
		product = powl(zz,2)*xx;
		break;
	case 19:
		product = powl(zz,2)*yy;
		break;
	case 20:
		product = xx*yy*zz;
		break;
	case 21:
		product = powl(xx,4);
		break;
	case 22:
		product = powl(yy,4);
		break;
	case 23:
		product = powl(zz,4);
		break;
	case 24:
		product = powl(xx,3)*yy;
		break;
	case 25:
		product = powl(xx,3)*zz;
		break;
	case 26:
		product = powl(yy,3)*xx;
		break;
	case 27:
		product = powl(yy,3)*zz;
		break;
	case 28:
		product = powl(zz,3)*xx;
		break;
	case 29:
		product = powl(zz,3)*yy;
		break;
	case 30:
		product = powl(xx,2)*pow(yy,2);
		break;
	case 31:
		product = powl(xx,2)*pow(zz,2);
		break;
	case 32:
		product = powl(yy,2)*pow(zz,2);
		break;
	case 33:
		product = powl(xx,2)*yy*zz;
		break;
	case 34:
		product = xx*powl(yy,2)*zz;
		break;
	case 35:
		product = xx*yy*powl(zz,2);
		break;
	default:
		break;
	}
	if (isnan(product)){
		std::cout <<"ERROR: NX:"<<NX << " xx: "<<xx <<" yy:" << yy << " zz: " << zz << " product: " << product << std::endl;
		product = 0.0;
		
	}
	ri = sqrt (pow(xx,2)+pow(yy,2)+pow(zz,2));
	//std::cout << " ri: " << ri << std::endl;
	double value = -ALPHA*pow(ri,2);
	//std::cout << " value: " << value << std::endl;
	double expValue = exp(value);
	//std::cout <<"xx: "<<xx <<" yy:" << yy << " zz: " << zz << " ALPHA: " << ALPHA << std::endl; 
	//std::cout << " expValue: " << expValue << std::endl;
	//std::cout << " pro: " << pro << std::endl;
	result = product * expValue;
	//std::cout << " Result: " << result << std::endl;
	if (isnan(result)){ 
		std::cout <<"ERROR: NX:"<<NX << " xx: "<<xx <<" yy:" << yy << " zz: " << zz << " result: " << result << std::endl;
		result = 0.0;
	}
	return result;
	
}
/**
 *  orb is a function to calculate the molecular orbital according to evaluation point
 * @param II
 * @param xx
 * @param yy
 * @param zz
 * @param NPO
 * @param nucleiData
 * @param IX
 * @param NX
 * @param ALPHA
 * @param C
 * @return 
 */
double orb (int II, double xx, double yy, double zz, int NPO, 
	std::vector <std::vector <double>> nucleiData, 
	std::vector <int>  IX, std::vector <int> NX,  std::vector <double> ALPHA, 
	std::vector <std::vector <double>> C){
	
	double orbResult=0.0;

	for (int i=0;i<NPO;i++){
		int NXi = NX.at(i);
		int IXi = IX.at(i) -1; //-1 vector element start in 0
		double X0i = nucleiData.at(IXi).at(1);
		double Y0i = nucleiData.at(IXi).at(2);
		double Z0i = nucleiData.at(IXi).at(3);
		double ALPHAi = ALPHA.at(i);
		double Ci = C.at(II).at(i);
		
		long double gauss = gaussian (NXi,xx-X0i,yy-Y0i,zz-Z0i,ALPHAi);
		/*std::cout << "ORB: xx: " << xx << " yy: " << yy << " zz: " << zz << " NXi: " << NXi << " IXi: " << IXi  
			<< " X0i: " << X0i  << " Y0i: " << Y0i << " Z0i: " << Z0i   
			<< " ALPHAi: " << ALPHAi  << " Ci: " << Ci  << " Gauss: " << gauss << std::endl;
		*/
		 orbResult = orbResult + Ci * gauss;	
	}
	if (isnan(orbResult)){ 
		std::cout << "ERROR orbResult: " << orbResult << std::endl;
		orbResult = 0.0;
	}
	return orbResult;	
}
/**
 * rho is a marginal density function
 * @param xx
 * @param yy
 * @param zz
 * @param NMO
 * @param NOCC
 * @param NPO
 * @param nucleiData
 * @param IX
 * @param NX
 * @param ALPHA
 * @param C
 * @return 
 */
long double rho (double xx, double yy, double zz, int NMO, std::vector <double> NOCC, int NPO, 
	std::vector <std::vector <double>> nucleiData, std::vector <int> IX, std::vector <int> NX, 
	std::vector <double> ALPHA, std::vector <std::vector <double>> C ){
	
	long double rhoResult = 0.0;
	long double orbi = 0.0;
	
	//std::cout << "RHO: xx: " << xx << " yy: " << yy << " zz: " << NMO  << " NPO: " << NPO  << std::endl;
	for (int i=0;i<NMO;i++){
		double NOCCi = NOCC.at(i);
		
		if (NOCCi> 0.0){
			orbi = orb (i,xx,yy,zz,NPO,nucleiData,IX,NX,ALPHA,C);
			rhoResult = rhoResult + NOCCi * pow(orbi,2);
		}
	}
	//std::cout << "rhoResult: " << rhoResult << std::endl;
	if (isnan(rhoResult)){
		std::cout << "ERROR: rhoResult: " << rhoResult << std::endl;
		rhoResult = 0.0;
	}
	return rhoResult;
	
}

/**
 * rho2 is a joint density of the evaluated point 
 * @param x1
 * @param y1
 * @param z1
 * @param x2
 * @param y2
 * @param z2
 * @param NMO
 * @param NOCC
 * @param numberA
 * @param NPO
 * @param nucleiData
 * @param IX
 * @param NX
 * @param ALPHA
 * @param C
 * @return 
 */
double rho2 (double x1, double y1, double z1, double x2, double y2, double z2, 
		int NMO, std::vector<double> NOCC, std::vector<double> numberA, 
		int NPO, std::vector<std::vector <double>> nucleiData, 
		std::vector<int> IX,std::vector<int> NX, std::vector<double> ALPHA, 
		std::vector <std::vector <double>> C ){
	
	double rho2Result = 0.0;
	
	double orbi1 = 0.0, orbi2 = 0.0, orbj1 = 0.0, orbj2 = 0.0;
	for (int i=0;i< NMO;i++){
		int NOCCi = NOCC.at(i);
		
		if (NOCCi>0.0){
			orbi1 = orb (i,x1,y1,z1,NPO,nucleiData,IX,NX,ALPHA,C);
			orbi2 = orb (i,x2,y2,z2,NPO,nucleiData,IX,NX,ALPHA,C);
			for (int j=0;j<NMO;j++){
				int NOCCj = NOCC.at(j);
				if (NOCCj>0.0){
					orbj1 = orb (j,x1,y1,z1,NPO,nucleiData,IX,NX,ALPHA,C);
					orbj2 = orb (j,x2,y2,z2,NPO,nucleiData,IX,NX,ALPHA,C);
					int numberAi = numberA.at(i);
					int numberAj = numberA.at(j);
					rho2Result = rho2Result + (NOCCi*NOCCj-numberAi*numberAj) * 
						(orbi1*orbj2*orbi1*orbj2 - orbi1*orbj2*orbj1*orbi2);
				}		
			}
		}
	}
	//std::cout << "rho2Result: " << rho2Result << std::endl;
	if (isnan(rho2Result)){
		std::cout << "ERROR: rho2Result: " << rho2Result << std::endl;
		rho2Result = 0.0;
	}
	return rho2Result;
}
/**
 * rho_cond is the conditional probability 
 * @param x1
 * @param y1
 * @param z1
 * @param x2
 * @param y2
 * @param z2
 * @param NMO
 * @param NOCC
 * @param numberA
 * @param NPO
 * @param nucleiData
 * @param IX
 * @param NX
 * @param ALPHA
 * @param C
 * @return 
 */
long double rho_cond  (double x1, double y1, double z1, double x2, double y2, double z2, 
		int NMO, std::vector<double>  NOCC, std::vector<double> numberA, int NPO, 
		std::vector<std::vector <double>> nucleiData, 
		std::vector<int>IX, std::vector<int> NX, std::vector<double> ALPHA, 
		std::vector <std::vector <double>> C ){
	
	long double rho_condResult = 0.0;
	
	rho_condResult = rho2(x1,y1,z1,x2,y2,z2,NMO,NOCC,numberA,NPO,nucleiData,IX,NX,ALPHA,C)/rho(x2,y2,z2,NMO,NOCC,NPO,nucleiData,IX,NX,ALPHA,C);
	//std::cout << "rho_condResult: " << rho_condResult << std::endl;
	if (isnan(rho_condResult)){
		std::cout << "ERROR: rho_condResult: " << rho_condResult << std::endl;
		rho_condResult = 0.0;
	}
	return rho_condResult;
}
/**
 * f_cutoff is a cutting function to restrict the range of calculation
 * @param rhoValue
 * @param rho_cut
 * @return 
 */
double f_cutoff (double rhoValue, double rho_cut){
	double f_cutoffResult = 0.0;
	
	f_cutoffResult = 0.5e0*(1.0e0 + erf(0.5e0 * log(rhoValue/rho_cut)));
	
	if (isnan(f_cutoffResult)){
		std::cout << "ERROR: f_cutoffResult: " << f_cutoffResult << std::endl;
		f_cutoffResult = 0.0;
	}
	return f_cutoffResult;
}


/**
 * KLD1_Func is a general function to calculate the values
 * @param x1
 * @param y1
 * @param z1
 * @param x2
 * @param y2
 * @param z2
 * @param NMO
 * @param NEL
 * @param NOCC
 * @param numberA
 * @param NPO
 * @param nucleiData
 * @param IX
 * @param NX
 * @param ALPHA
 * @param C
 * @param rho_cut
 * @return 
 */
long double KLD1_Func (double x1, double y1, double z1, double x2, double y2, double z2,
		int NMO, double NEL, std::vector <double> NOCC, std::vector <double> numberA, int NPO, 
		std::vector <std::vector <double>> nucleiData, std::vector <int> IX, 
		std::vector <int> NX, std::vector <double> ALPHA,
		std::vector <std::vector <double>> C, double rho_cut){
	
	//std::cout << "X1: " << x1 << " Y1: " << y1 << " Z1: "  << z1 << " X2: " << x2 << " Y2: " << y2 << " Z2: "  << z2 << std::endl;
	
	
	long double resultKLD1_Func = 0.0;
	//std::cout << "Starting calculation..." << std::endl;
	long double rcd = rho_cond(x1,y1,z1,x2,y2,z2,NMO,NOCC,numberA,NPO,nucleiData,IX,NX,ALPHA,C)/(NEL-1.e0);
	//std::cout << "rcd: "<<rcd<<std::endl;
	double rx1 = rho(x1,y1,z1,NMO,NOCC,NPO,nucleiData,IX,NX,ALPHA,C)/NEL;
	//std::cout << "rx1: "<<rx1<<std::endl;
	double fco = f_cutoff(rho(x2,y2,z2,NMO,NOCC,NPO,nucleiData,IX,NX,ALPHA,C),rho_cut);
	//std::cout << "fco: "<<fco<<std::endl;
	
	if ((rcd <= 0.0)|| (rx1<=0.0)){
		resultKLD1_Func = 0.0;
	}else{
		resultKLD1_Func = (NEL-1.e0 ) * fco * rcd * (log(rcd/rx1)/log(2.e0));
	}
	if (isnan(resultKLD1_Func)){
		std::cout << "ERROR: resultKLD1_Func: " << resultKLD1_Func << std::endl;
		resultKLD1_Func = 0.0;
	}
	
	//std::cout << "KLD1D= " << resultKLD1_Func << std::endl;
	return resultKLD1_Func;
	
}



/**
 * KLD2_Func is a calculation if the number of wfn is two
 * @param x1
 * @param y1
 * @param z1
 * @param x2
 * @param y2
 * @param z2
 * @param NMO1
 * @param NEL1
 * @param NOCC1
 * @param NA1
 * @param NPO1
 * @param X01
 * @param Y01
 * @param Z01
 * @param IX1
 * @param NX1
 * @param ALPHA1
 * @param C1
 * @param NMO2
 * @param NOCC2
 * @param NEL2
 * @param NA2
 * @param NPO2
 * @param X02
 * @param Y02
 * @param Z02
 * @param IX2
 * @param NX2
 * @param ALPHA2
 * @param C2
 * @param rho_cut
 * @return 
 */
double KLD2_Func (double x1, double y1, double z1, double x2, double y2, double z2, 
	int NMO1, double NEL1,std::vector<double> NOCC1, std::vector<double> numberA1, 
	int NPO1, std::vector<std::vector<double>> nucleiData1, std::vector<int> IX1, 
	std::vector<int> NX1,std::vector<double> ALPHA1, std::vector <std::vector <double>> C1,
	int NMO2, std::vector<double> NOCC2, double NEL2, std::vector<double> numberA2, 
	int NPO2, std::vector<std::vector<double>> nucleiData2, std::vector<int> IX2, 
	std::vector<int> NX2, std::vector<double> ALPHA2, std::vector <std::vector <double>> C2,  
	double rho_cut){
	
	double resultKLD2_Func = 0.0;
	
	double rcd1 = rho_cond(x1,y1,z1,x2,y2,z2,NMO1,NOCC1,numberA1,NPO1,nucleiData1,IX1,NX1,ALPHA1,C1)/(NEL1-1.e0);
	double rcd2 = rho_cond(x1,y1,z1,x2,y2,z2,NMO2,NOCC2,numberA2,NPO2,nucleiData2,IX2,NX2,ALPHA2,C2)/(NEL2-1.e0);
	double fco = f_cutoff(rho(x2,y2,z2,NMO1,NOCC1,NPO1,nucleiData1,IX1,NX1,ALPHA1,C1),rho_cut);
	if ((rcd1 <= 0.0)|| (rcd2<=0.0)){
		resultKLD2_Func = 0.0;
	}else{
		resultKLD2_Func = NEL1 * fco * rcd1 * (log(rcd1/rcd2)/log(2.0e0));
	}
	
	return resultKLD2_Func;
}
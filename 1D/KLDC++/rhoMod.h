/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   rhoMod.h
 * Author: vrodriguez
 *
 * Created on September 24, 2018, 3:57 PM
 */
#include <vector>

#ifndef RHOMOD_H
#define RHOMOD_H
double gaussian (int NX, double xx, double zz, double ALPHA);
double orb (int II, double xx, double yy, double zz, int NPO, 
	std::vector <std::vector <double>> nucleiData, 
	std::vector <int>  IX, std::vector <int> NX,  std::vector <double> ALPHA, 
	std::vector <std::vector <double>> C);
long double rho (double xx, double yy, double zz, int NMO, std::vector <double> NOCC, int NPO, 
	std::vector <std::vector <double>> nucleiData, std::vector <int> IX, std::vector <int> NX, 
	std::vector <double> ALPHA, std::vector <std::vector <double>> C );
double rho2 (double x1, double y1, double z1, double x2, double y2, double z2, 
		int NMO, std::vector<double> NOCC, std::vector<double> numberA, 
		int NPO, std::vector<std::vector <double>> nucleiData, 
		std::vector<int> IX,std::vector<int> NX, std::vector<double> ALPHA, 
		std::vector <std::vector <double>> C );
long double rho_cond  (double x1, double y1, double z1, double x2, double y2, double z2, 
		int NMO, std::vector<double>  NOCC, std::vector<double> numberA, int NPO, 
		std::vector<std::vector <double>> nucleiData, 
		std::vector<int>IX, std::vector<int> NX, std::vector<double> ALPHA, 
		std::vector <std::vector <double>> C );
double f_cutoff (double rhoValue, double rho_cut);
long double KLD1_Func (double x1, double y1, double z1, double x2, double y2, double z2,
		int NMO, double NEL, std::vector <double> NOCC, std::vector <double> numberA, int NPO, 
		std::vector <std::vector <double>> nucleiData, std::vector <int> IX, 
		std::vector <int> NX, std::vector <double> ALPHA,
		std::vector <std::vector <double>> C, double rho_cut);
#endif /* RHOMOD_H */

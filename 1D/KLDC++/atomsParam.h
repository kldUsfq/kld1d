/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   atomsParam.h
 * Author: vrodriguez
 *
 * Created on September 17, 2018, 10:01 PM
 */
#include <string>


 /**
  * returns the symbol according to the number
  * @param symbolNumber
  * @return 
  */
 std::string getSymbol (int symbolNumber);
/**
 * returns the nbc for symbolNumber
 * @param symbolNumber of the symbol
 * @return 
 */
int getNbc (int symbolNumber);
/**
 * returns the  bsr of symbolNumber
 * @param symbolNumber
 * @return 
 */
double getBsr (int symbolNumber);
/**
 * finds the symbol and returns the position
 * @param symbol to search
 * @return position
 */
int getSymbolNumber (std::string symbol);






/* 
 * File:   common.h
 * Author: vrodriguez
 *
 * Created on September 15, 2018, 9:03 PM
 */

#ifndef COMMON_H
#define COMMON_H

#include <string>
//
//  common.h
//  KLDC++
//
//  Created by Vladimir Rodriguez on 8/27/18.
//  Copyright © 2018 Universidad San Francisco de Quito. All rights reserved.
//
using namespace std;
bool fileExists (const string& name);
int paramsValidation (const char * app[], const char  * kldV[], const char * paramFile[]);
int wfn_input (string paramFile);


#endif /* COMMON_H */


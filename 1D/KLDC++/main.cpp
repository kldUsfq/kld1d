//
//  main.cpp
//  KLD_Intel_CPP
//  Version : 1.0.0  (Fortran 90/MPI/OpenMP)
//  Author: Luis Rincon, USFQ, Aug 27, 2015
//  Revision: Luis Rincon, USFQ, Sep 10,2017
//  Version : 1.0.0 (c++/MPI/OpenMP)
//  Vladimir Rodriguez on 8/27/18.
//  Copyright © 2018 Universidad San Francisco de Quito. All rights reserved.
//

#include <iostream>
#include <time.h>
#include <ctime>
#include <chrono>
#include "sharedkld.h"

using namespace std;

/**
 * Variables
 */
//path where the application is run
string appPath;
//version of KLD: 1D, 2D, 3D
string kldVersion;
//file that contains the parameters for running KLD
string paramsPath;

/**
 * KLD 1D method
 * calculate the Kullback-Leibler divergence in a 1D-grid
 * Version : 1.0.0  (Fortran 90/MPI/OpenMP)
 * Author: Luis Rincon, USFQ, Aug 27, 2015
 * Revision: Luis Rincon, USFQ, Sep 10,2017
 * C++ Author: Vladimir Rodriguez, USFQ, Sep 2018
 */
int KLD_1D (string parameters){
    //time counting  variables
    float  cpu_t1=0,cpu_t2=0;
    using chrono::system_clock;
    std::chrono::system_clock::time_point startTime;
    std::chrono::system_clock::time_point endTime;
    
    //get cpu time at start
    cpu_t1 = time(NULL);
    
    cout << "Initializing calculations for KLD_1D..." << endl;
    // Start time
    startTime = chrono::system_clock::now();
    std::time_t timeS;
    timeS = std::chrono::system_clock::to_time_t(startTime);
    cout << "Time at start is: " << ctime (&timeS) << endl;
    //staring calculations
    //read input parameters from specified file
    int result1D = wfn_input (parameters);
    
    //get cpu time at end
    cpu_t2 = time (NULL);
    endTime = chrono::system_clock::now();
    std::time_t timeE;
    timeE = std::chrono::system_clock::to_time_t(endTime);
    cout << "Time at end is: " << ctime(&timeE) << endl;
    cout << "Elapsed time is: " << std::chrono::system_clock::to_time_t(endTime) - std::chrono::system_clock::to_time_t(startTime) << endl;
    cout << "Runtime: " << cpu_t2 - cpu_t1 << " seconds." << endl;
    return 0;
}

/**
 * Main Function of KLD
 */
int main(int argc, const char * argv[]) {
    
    //parsing arguments from command line
    if (argc == 3){
        if (paramsValidation (&argv [0], &argv[1], &argv[2])==0){
            appPath = argv [0];
            kldVersion = argv[1];
            paramsPath = argv[2];
            cout << "Using the following arguments:" << endl;
            cout << "Path of APP:" << appPath << endl;
            cout << "KLD version:" << kldVersion << endl;
            cout << "Parameters file:" << paramsPath << endl;
            /**
             * Calculating 1D
             */
            if (kldVersion=="1D"){
                cout << "Processing 1D version..." << endl;
                int success = KLD_1D (paramsPath);
                if (success == 0){
                    cout << "End processing 1D."<< endl;
                }else {
                    cout << "Fail processing 1D."<< endl;
                }
            /**
            * Calculating 2D
            */
            }else if (kldVersion=="2D"){
                cout << "Processing 2D version..." << endl;
            /**
             * Calculating 3D
            */
            }else if (kldVersion=="3D"){
                cout << "Processing 3D version..." << endl;
            }
            
            
        } else {
            cout << "Not valid arguments specified: 1D,2D, 3D for kld version" << endl << "Path of the file where parameters are specified: i.e /home/user/file.in" << endl;
        }
    }else{
        cout << "Not valid arguments specified: 1D,2D, 3D for kld version" << endl << "Path of the file where parameters are specified: i.e /home/user/file.in" << endl;
     }
    return 0;
}




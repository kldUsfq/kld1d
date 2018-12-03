//
//  common.cpp
//  KLD_Intel_CPP
//  this module read the wfn input file from Gaussian or Gamess
//  Author: Luis Rincon, ULA, 03.2014 (Include wfn_type :: 24/07/2017)

//  C++ version created by Vladimir Rodriguez on 8/27/18.
//  Copyright © 2018 Universidad San Francisco de Quito. All rights reserved.
/*
*@TODO: DEFINE A STRUCTURE/CLASS  WITH VECTOR FIELDS INSTEAD OF INDIVIDUAL FIELDS
WFN INPUT SHOULD LOAD THE DATA OF CLASS WFN
*/


#include <iostream>
#include <stdio.h>
#include <fstream>
#include <string>
#include <dirent.h>
#include <regex>
#include <math.h>
#include <sstream>
//#include <Accelerate/Accelerate.h>
#include "/opt/intel/mkl/include/mkl_lapack.h"
#include <numeric>
#include "sphere_lebedev_rule.h"
#include "atomsParam.h"
#include "rhoMod.h"
#include "/opt/intel/include/omp.h"

using namespace std;

/**
 * # of WFN that will be process
 */
string numberWFN;

/**
 * variables to determine size of the matrix
 */
int rowsMol = 0;
int columsPrim = 0;
int nuclei = 0;
string method = "";
/**
 number of electrons
 vector to determine the error
 */
vector <double> numberA;
double numberElectrons=0.0;
/**
 *vector for nuclear center for each primitive
 */
vector <int> nucCenterData;

/**
 *vector for orbital  for each primitive
 */
vector <int> orbTypeData;

/**
 *vector for orbital exponent for each primitive
 */
vector <double> orbExpData;

/** vector for ocuppational number for each molecular orbital*/
vector <double> ocupNum;

/**
 matrix 2d for the data related to nuclei
 */
vector <vector <double>> nucleiData;
vector <double> rowNucleiData;

/**
 * matrix 2D for the data with dimmensions of the array Molecular Orbital
 */
vector <vector <double>> molOrbiData;
vector <double> rowMolOrbiData;

/**
 * prototype for debug and print results of DGESV
 */
extern void print_matrix( string desc, int m, int n, double* a, int lda );
extern void print_int_vector( string desc, int n, int* a );

//name of the file where results will be saved
string resultsFile;


/**
 * Levenberg Marqueardt parameters
 */

double epsilon1 = 1e-3;
double epsilon2 = 1e-3;
double epsilon3 = 1e-7;
double epsilon4 = 1e-1;
double lambda0 = 1e-2;
double lambda_up = 11.0;
double lambda_dn = 9.0;


/**
 * funtion that executes Levenberg Marquardt optimization of chiˆ2
*/
int optLevMarOpt (){
    //variables for DGESV Calculation
    int n = ocupNum.size();
    int nrhs = 1;
    int info;
    int LDA = ocupNum.size();
    
    int LDB = ocupNum.size();
    int IPIV[n];
    
    //max number of iterations
    int maxIter = (int) 10*ocupNum.size();
    
    vector <double> YI;
    vector <double> YJ;

    //initialization of the parameters
    double chi2n = 0.0;
    for (const double& s1 : ocupNum) {
        double yi = s1 * (s1 - 1.0);
        YI.push_back(yi);
        YJ.push_back(0.0);
    }
    for (int i=0;i< ocupNum.size();i++) {
        for (int j=0;j<ocupNum.size();j++) {
            if (i != j){
                YJ.at(i)=YJ.at(i)-numberA.at(i)*numberA.at(j);
            }
        }
        chi2n = chi2n + pow(YI.at(i)-YJ.at(i),2);
        
    }
    double chi2o = 0.0;
    double lambda = lambda0;
    
    bool iconv = true;
    if (chi2n < epsilon3){
        iconv = false;
    }
    int niter = 0;
    
    // processing and optimization
    while (iconv){
        niter++;
        cout << "LM Iteration: "<< niter << " CHI2: "<< chi2n << " DIF: " << chi2n - chi2o << endl;
        chi2o = chi2n;
	
        /**
         * Jacobian vector 2D and initializacion
         */
	/**
	* @TODO make th is part of jacobian parallelized 
	*/
        vector <vector <double>> jacobianData;
        vector <double> jacobianRow;
        for (int i=0;i<ocupNum.size();i++){
            for (int j=0;j<ocupNum.size();j++){
                jacobianRow.push_back(0.0);
            }
            jacobianData.push_back(jacobianRow);
            jacobianRow.clear();
        }
	
	
        for (int i=0;i<ocupNum.size();i++){
            for (int j=0;j<ocupNum.size();j++){
                if (i==j){
                    jacobianData.at(i).at(i) = 0;
                    for (int k=0;k<ocupNum.size();k++){
                        if (i!=k){

                            jacobianData.at(i).at(i) = jacobianData.at(i).at(i) - numberA.at(k);
                        }
                    }
                }else {
                    jacobianData.at(i).at(j) = -1 * numberA.at(i);
                }
            }
        }
	
        /**
         * Hessian
         */
         
        //initialization of hessian vector 2D
        vector <vector <double>> hessianData;
        vector <double> hessianRow;
        for (int i=0;i<ocupNum.size();i++){
            for (int j=0;j<ocupNum.size();j++){
                hessianRow.push_back(0.0);
            }
            hessianData.push_back(hessianRow);
            hessianRow.clear();
        }
        //data procesing for hessian
	//#pragma omp for
	
        for (int i=0;i<ocupNum.size();i++){
            for (int j=0;j<ocupNum.size();j++){
                hessianData.at(i).at(j) = 0;
                for (int k=0;k<ocupNum.size();k++){
                    hessianData.at(i).at(j) = hessianData.at(i).at(j) + jacobianData.at(k).at(i)*jacobianData.at(k).at(j);
                }
            }
        }
	
        //initialization of A and B 2D vectors
        vector <vector <double>> AData;
        vector <double> ARow;
        vector <vector <double>> BData;
        vector <double> BRow;
        for (int i=0;i<ocupNum.size();i++){
            for (int j=0;j<ocupNum.size();j++){
                ARow.push_back(0.0);
                BRow.push_back(0.0);
            }
            AData.push_back(ARow);
            ARow.clear();
            BData.push_back(BRow);
            BRow.clear();
        }
        //processing for hessian A 2D vector
	//#pragma omp for
        for (int i=0;i<ocupNum.size();i++){
            for (int j=0;j<ocupNum.size();j++){
                if (i==j){
                    AData.at(i).at(i)=(1.0+lambda)*hessianData.at(i).at(i);
                }else {
                    AData.at(i).at(j)=hessianData.at(i).at(j);
                    
                }
                
            }
        }
        //processing for hessian B 2D vector
        for (int i=0;i<ocupNum.size();i++){
            for (int j=0;j<ocupNum.size();j++){
                BData.at(i).at(0)=BData.at(i).at(0)+jacobianData.at(j).at(i)*(YI.at(j)-YJ.at(j));
            }
        }
        //copy the elements from the 2D vector A to linear array a
        double a[LDA*n];
        int elementA = 0;
        for (int i=0;i<LDA;i++){
            for (int j=0;j<n;j++){
                a[elementA]=AData.at(i).at(j);
                elementA++;
            }
        }
        int elementB = 0;
        double b[LDB*nrhs];
        //copy the elements from the 2D vector B to linear array b
        for (int i=0;i<LDB;i++){
            for (int j=0;j<nrhs;j++){
                b[elementB]=BData.at(i).at(j);
                elementB++;
            }
        }
        /* Computes the solution to the system of linear
        * equations with a square matrix A and multiple
        * right-hand sides B, where A is the coefficient matrix
        * and B is the right-hand side matrix. USES INTEL LAPACK
        */
        
        dgesv_(&n, &nrhs,a, &LDA, IPIV, b, &LDB, &info);
        if (info > 0){
            cout << "DGESV ERROR: Cannot find a solution"<< endl;
            break;
        }
        /*
         FOR  DEBUG
         */
        //print_matrix("Solution:", n, nrhs, b, LDB);
        //print_matrix("Details of LU factorization: ", n,n, a, LDA);
        //print_int_vector("Pivot indices:", n, IPIV);
        
        //convert a 1D to A 2D
        elementA = 0;
        for (int j=0;j<LDA;j++){
            for (int k=0;k<n;k++){
                AData.at(j).at(k) = a[elementA];
                elementA++;
            }
        }
        //convert b 1D to B 2D
        elementB = 0;
        for (int j=0;j<LDB;j++){
            for (int k=0;k<nrhs;k++){
                BData.at(j).at(k) = b[elementB];
                elementB++;
            }
        }
        // new chi2n
        chi2n = 0.0;
        double yyj = 0.0;
        for (int i=0;i<ocupNum.size();i++){
            yyj = 0.0;
            for (int j=0;j<ocupNum.size();j++){
                if (i!=j){
                    yyj = yyj - (numberA.at(i) + BData.at(i).at(0)) * ((numberA.at(j) + BData.at(j).at(0)));
                }
            }
            chi2n = chi2n + pow(YI.at(i)-yyj,2);
        }
        //RHOI
        double rhoi = 0.0;
        double den1 = 0.0;
	#pragma omp for
        for (int i=0;i<ocupNum.size();i++){
            den1 = den1 + lambda*BData.at(i).at(0) * hessianData.at(i).at(i) * BData.at(i).at(0);
        }
        double den2 = 0.0;
	
        for (int i=0;i<ocupNum.size();i++){
	    #pragma omp for
            for (int j=0;j<ocupNum.size();j++){
                den2 = den2 + BData.at(i).at(0) * jacobianData.at(j).at(i) * (YI.at(j)-YJ.at(j));
            }
        }
        rhoi = (chi2o-chi2n)/(den1+den2);
        
        //ACEPTANCE CRITERIA
        if (rhoi>epsilon4){
            lambda = std::max(lambda/lambda_dn,1e-7);
            for (int i=0;i<ocupNum.size();i++){
                numberA.at(i)= numberA.at(i) + BData.at(i).at(0);
            }
            for (int i=0;i<ocupNum.size();i++){
                YJ.at(i)=0.0;
                for (int j=0;j<ocupNum.size();j++){
                    if (i!=j){
                        YJ.at(i)= YJ.at(i) - numberA.at(i)*numberA.at(j);
                    }
                }
            }
        }else{
            lambda=std::min(lambda*lambda_up,1e+7);
        }
        //CONVERGENCE CRITERIA
        //gradient
        double max_grad = 0.0;
        double grad = 0.0;
        for (int i=0;i<ocupNum.size();i++){
            grad=0.0;
            for (int j=0;j<ocupNum.size();j++){
                grad = grad + jacobianData.at(j).at(i)*(YI.at(j)-YJ.at(j));
            }
            if (fabs(grad)>max_grad){
                max_grad=fabs(grad);
            }
        }
        //PARAMETERS
        double max_par = 0.0;
        double par=0.0;
        for (int i=0;i<ocupNum.size();i++){
            par = BData.at(i).at(0)/YJ.at(i);
            if (fabs(par)>max_par){
                max_par = fabs(par);
            }
        }
        if (max_grad<epsilon1){
            iconv=false;
        }
        if (max_par<epsilon2){
            iconv=false;
        }
        if (chi2n<epsilon3){
            iconv=false;
        }
        if (niter>maxIter){
            iconv=false;
        }

    }//while (iconv)
    
    cout << "CHIˆ2 = " << chi2n << endl;
    return chi2n;
}


///variables for cube
double cubeX, cubeY,cubeZ;
double cubeXData [4];
double cubeYData [4];
double cubeZData [4];

/**
 * process cube.dat file and loads information from there.
 * @return 
 */
string cubeInput (string paramFile){
	string cubeOutputFileName = "";
	std::size_t found = paramFile.find_last_of("/\\");
	string cubeFilePath = paramFile.substr(0,found+1);
    
	cubeFilePath.append("cube.dat");
	cout << "Reading data file:  "<< cubeFilePath << "..." <<endl;
    	std::ifstream dataPath(cubeFilePath.c_str());
	int line =1;
	regex reg ("\\s+");
	std::string str;
	bool loadNewRow = false;
	while (std::getline(dataPath, str)){
		//tokenized for str
		sregex_token_iterator iter(str.begin(), str.end(), reg, -1);
		sregex_token_iterator end;
		vector<string> fieldsStr (iter, end);
		if (line==1){
			string resultFilePath = paramFile.substr(0,found+1);
			//removes ' ' from the name of the file
			str = str.substr(1,str.length()-2);
			resultFilePath.append(str);
			cubeOutputFileName = resultFilePath;
		}else if (line ==2){
			//reads coordinates x,y,z of first point in the cube
			std::istringstream os1(fieldsStr[1]);
			os1 >> cubeX; // is XX0
			std::istringstream os2(fieldsStr[2]);
			os2 >> cubeY; // is YY0
			std::istringstream os3(fieldsStr[3]);
			os3 >> cubeZ; // is ZZ0
		}else if (line ==3){
			//reads number of point in X
			std::istringstream os1(fieldsStr[0]);
			os1 >> cubeXData[0]; // is INX
			std::istringstream os2(fieldsStr[1]);
			os2 >> cubeXData[1]; // is IXA
			std::istringstream os3(fieldsStr[2]);
			os3 >> cubeXData[2]; // is IXB
			std::istringstream os4(fieldsStr[3]);
			os4 >> cubeXData[3]; // is IXC
		}else if (line ==4){
			//reads number of point in Y
			std::istringstream os1(fieldsStr[0]);
			os1 >> cubeYData[0];
			std::istringstream os2(fieldsStr[1]);
			os2 >> cubeYData[1];
			std::istringstream os3(fieldsStr[2]);
			os3 >> cubeYData[2];
			std::istringstream os4(fieldsStr[3]);
			os4 >> cubeYData[3];
		}else if (line ==5){
			//reads number of point in Z
			std::istringstream os1(fieldsStr[0]);
			os1 >> cubeZData[0];
			std::istringstream os2(fieldsStr[1]);
			os2 >> cubeZData[1];
			std::istringstream os3(fieldsStr[2]);
			os3 >> cubeZData[2];
			std::istringstream os4(fieldsStr[3]);
			os4 >> cubeZData[3];
		}
		line++;
	}
	dataPath.close();
	cout << "Read cube.dat OK ... " << line-1 << " lines" << endl; 
	cout << "Results will be saved in: " << cubeOutputFileName << " file" << endl;
	
	
	return cubeOutputFileName;
}

/**
 * WEIGHT OF IA ATOM AT THE POSTION X,Y,Z (GENERALIZED VORONOI POLIHEDRA) PROPOSED BY BECKE
 * THIS WEIGHT IS RANGE BETWEEN 0 AND 1
 * @param i
 * @param x
 * @param y
 * @param z
 * @return 
 */
double wv (int ia, double x, double y, double z){
	double resultWV = 0.0;
	
	const int mu_max = 3;
	double ri =0.0;
	double rj = 0.0;
	double Rij = 0.0;
	double muij = 0.0, chij = 0.0,uij = 0.0, aij = 0.0;
	double *ww = new double [nuclei];
	//initializes the array
	for (int k=0;k<nuclei;k++){
		ww[k]=1.0;
	}
	
	for (int i=1;i<nuclei;i++){
		for (int j=0;j<=i-1;j++){
			double nucleiDataX = nucleiData.at(i).at(1);
			double nucleiDataY = nucleiData.at(i).at(2);
			double nucleiDataZ = nucleiData.at(i).at(3);
			double nuclejDataX = nucleiData.at(j).at(1);
			double nuclejDataY = nucleiData.at(j).at(2);
			double nuclejDataZ = nucleiData.at(j).at(3);
			ri = sqrt ( pow(x - nucleiDataX,2) + pow(y - nucleiDataY,2) + pow(z - nucleiDataZ,2) ); 
			
			
			rj = sqrt ( pow(x - nuclejDataX,2) + pow(y - nuclejDataY,2) + pow(z - nuclejDataZ,2) ); 
			Rij = sqrt (pow(nucleiDataX-nuclejDataX,2) + pow(nucleiDataY-nuclejDataY,2) + pow(nucleiDataZ-nuclejDataZ,2));
			muij = (ri - rj)/Rij;
			
			//estimation of aij
			double rbsi = getBsr (nucleiData.at(i).at(0));
			
			if (getSymbol(nucleiData.at(i).at(0))!="H"){
				rbsi = 2.0 * rbsi;
			}
			
			double rbsj = getBsr (nucleiData.at(j).at(0));
			if (getSymbol(nucleiData.at(j).at(0))!="H"){
				rbsj = 2.0 * rbsj;
			}
			
			chij = rbsi/rbsj;
			uij = (chij-1.0)/(chij+1.0);
			aij = uij/(pow(uij,2)-1.0);
			
			if (aij > 0.5){
				aij = 0.5;
			}
			if (aij < -0.5){
				aij = -0.5;
			}
			muij = muij + aij * (1.0 - pow(muij,2));
				
			for (int k=0;k<mu_max;k++){
				muij = 1.5 * muij -0.5 * pow(muij,3);
			}
			
			ww[i]=0.5 * ww[i] * (1-muij);
			ww[j]=0.5 * ww[j] * (1+muij);
			
		}
		
		
	}
	double sumAcum = 0.0;
	
	for (int counter=0;counter<nuclei;counter++){
		sumAcum = sumAcum + ww[counter];
	}
	resultWV = ww[ia]/sumAcum;
	
	//release reserved memory
	delete [] ww;
	return resultWV;
}


//variables for becke
int rule_lebedev=7;
double *totalInt, *partial_int;
int available_lebedev, order_lebedev, precision_lebedev;
double  *x_lebedev,*y_lebedev,*z_lebedev,*w_lebedev;
long double *x,*y,*z,*w_becke;
int npoints,kpoint;
long double KLD = 0.0;

double X1= 0.0, Y1= 0.0, Z1 = 0.0;

double domeSigma = 4.5;
double r12 = 0.0;
const double rho_cut=0.001e0;
/**
 * NUMERICAL INTEGRATION PROGRAM (MPI/OPENMP)
 * @return si the process was ok
 */
bool becke88 (){
	bool calc = false;
	KLD = 0.0;
	
	/* call to the subroutines that generate Lebedev grids for integration on a sphere.
	!  the subroutines are included in the file SPHERE_LEBEDEV_RULE
	!
	!  Author:
	!
	!    Dmitri Laikov
	!
	!  Reference:
	!
	!    Vyacheslav Lebedev, Dmitri Laikov,
	!    A quadrature formula for the sphere of the 131st
	!    algebraic order of accuracy,
	!    Russian Academy of Sciences Doklady Mathematics,
	!    Volume 59, Number 3, 1999, pages 477-481.
	!
	*/
	available_lebedev = available_table (rule_lebedev);
	if (available_lebedev==1){
		order_lebedev = order_table (rule_lebedev);
		cout << "order lebedev: "<< order_lebedev << endl;
		precision_lebedev = precision_table (rule_lebedev);
		cout << "precision lebedev: "<< precision_lebedev << endl;
		//init the arrays
		x_lebedev = new  double[order_lebedev];
		y_lebedev = new  double[order_lebedev];
		z_lebedev = new  double[order_lebedev];
		w_lebedev = new  double[order_lebedev];
		//call method LD_BY_ORDER returns a Lebedev angular grid given its order.
		ld_by_order (order_lebedev, x_lebedev, y_lebedev, z_lebedev, w_lebedev);

	}
	npoints  = 0;
	for (int i=0;i<nuclei;i++){
		npoints = npoints + getNbc(nucleiData.at(i).at(0)) * order_lebedev;
	}
	cout << "Number of points: " << npoints << endl;
	
	//initializes x,y,z and w_becke arrays with the number of points
	x = new long double[npoints];
	y = new long double[npoints];
	z = new long double[npoints];
	w_becke = new long double[npoints];
	int kpoint = -1;
	for (int i=0;i<nuclei;i++){
	    //Bragg-Slater radius of atom ia
		double r_bs = getBsr (nucleiData.at(i).at(0));
		// radial loop: Gauss-Chevichev
		double nbc = getNbc(nucleiData.at(i).at(0));
		int j =0, k=0;
		
		for (j=1;j<=nbc;j++){
			//Gauss-Chevichev quadrature grid
			long double x_gc = cos ((j/(nbc+1)) * M_PI) ;
			long double w_gc = (M_PI / (nbc+1)) * pow(sin((j/(nbc+1)) * M_PI),2);
			long double wgc = 2.0 * pow(r_bs,3) * w_gc * sqrt (pow(1.0+x_gc,3) / pow(1.0-x_gc,9));

			// Gauss-Chevichev radius
			long double r_gc = r_bs * ((1.0+x_gc)/(1.0-x_gc));
			//angular loop: Levedevs quadrature
			for (k=0;k<order_lebedev;k++){
				kpoint = kpoint+1;
				double xi = nucleiData.at(i).at(1);
				double yi = nucleiData.at(i).at(2);
				double zi = nucleiData.at(i).at(3);
				//nucleiData.at(i).at(1) is nucleiData [i,1] where i is row and 1 is column that refer to X
				x[kpoint] =  xi + x_lebedev[k] * r_gc;
				//nucleiData.at(i).at(1) is nucleiData [i,2] where i is row and 1 is column that refer to Y
				y[kpoint] =  yi + y_lebedev[k] * r_gc;
				//nucleiData.at(i).at(1) is nucleiData [i,3] where i is row and 1 is column that refer to Z
				z[kpoint] =  zi + z_lebedev[k] * r_gc;
				w_becke[kpoint] = 4.0 * M_PI * w_lebedev[k] * wv(i,x[kpoint],y[kpoint],z[kpoint]) * wgc;
				
			}
		}
	}
	
	/*for (int cont=0;cont<100;cont++){
		cout << "X[" << cont << "]: " << x[cont];
		cout << " Y[" << cont << "]: " << y[cont];
		cout << " Z[" << cont << "]: " << z[cont] << endl;
		
	}
	exit (0);
	*/
	//releases reserved memory
	delete [] x_lebedev;
	delete [] y_lebedev;
	delete [] z_lebedev;
	delete [] w_lebedev;
	
	//initializes the vector partial_int of data according to conditions below
	if (cubeXData[0]==0 && cubeYData[0]==0){
		double size = cubeZData[0];
		partial_int = new double[size];
		for (int i=0;i<cubeZData[0];i++){
			partial_int[i]=0.0;
		}
	}else if (cubeXData[0]==0 && cubeZData[0]==0){
		double size = cubeYData[0];
		partial_int = new double[size];
		for (int i=0;i<cubeYData[0];i++){
			partial_int[i]=0.0;
	 	}
	}else if (cubeYData[0]==0 && cubeZData[0]==0){
		double size = cubeXData[0];
		partial_int = new double[size];
		for (int i=0;i<cubeXData[0];i++){
			partial_int[i]=0.0;
		}
	}
	int IN1 = 0;
	if ((cubeXData[0]==0)&&(cubeYData[0]==0)){
		IN1 = cubeZData[0];
	} 
	if ((cubeXData[0]==0)&&(cubeZData[0]==0)){
		IN1 = cubeYData[0];
	} 
	if ((cubeYData[0]==0)&&(cubeZData[0]==0)){
		IN1 = cubeXData[0];
	}
	
	/*
	 * @TODO proceso paralelo de cálculo y distribuido
	 */
	int nThreads = 0;
	const int factorN = 20;
	int sysThreadLimit = 0;
	
	//initializes the vector total of data according to conditions below
	
	
	if (cubeXData[0]==0 && cubeYData[0]==0){
		double size = cubeZData[0];
		totalInt = new  double[size];
		for (int i=0;i<cubeZData[0];i++){
			totalInt[i]=0.0;
		}
	}else if (cubeXData[0]==0 && cubeZData[0]==0){
		double size = cubeYData[0];
		totalInt = new  double[size];
		for (int i=0;i<cubeYData[0];i++){
			totalInt[i]=0.0;
		}
	}else if (cubeYData[0]==0 && cubeZData[0]==0){
		double size = cubeXData[0];
		totalInt = new  double[size];
		for (int i=0;i<cubeXData[0];i++){
			totalInt[i]=0.0;
		}
	}	
	
	
	
	#pragma omp parallel
	{	
		nThreads = omp_get_num_threads();
		sysThreadLimit = omp_get_thread_limit();
		if (omp_get_thread_num() == 0){
			std::cout << "Number physical cores of node:" << nThreads << std::endl;
			std::cout << "Number of Threads:" << nThreads*factorN << std::endl;
			std::cout << "System thread limit :" << sysThreadLimit << std::endl;
		}
			
	}
	int progress = 1;
	int IJK=0, III =0; 
	for ( IJK=0;IJK<npoints;IJK++){
		#pragma omp parallel for private (X1,Y1,Z1,III, KLD) num_threads (nThreads*factorN)
		for ( III =0;III<IN1;III++){
			
			
			if (cubeXData[0]==0){
				X1 = cubeX;
			} else {
				X1 = cubeX + III*cubeXData[1];
			}
			if (cubeYData[0]==0){
				Y1 = cubeY;
			}else{
				Y1 = cubeY + III*cubeYData[2];
			}
			if (cubeZData[0]==0){
				Z1 = cubeZ;
			}else{
				Z1 = cubeZ + III*cubeZData[3];
			}
			r12 = sqrt (pow((x[IJK]-X1),2) + pow((y[IJK]-Y1),2) + pow((z[IJK]-Z1),2));
			//excludes distant elements
			if (r12 <= domeSigma){
				if (numberWFN == "1"){
					KLD = KLD1_Func (x[IJK],y[IJK],z[IJK],
					         X1,Y1,Z1,
						 ocupNum.size(),numberElectrons, ocupNum,
						 numberA, orbTypeData.size(),
						 nucleiData,
						 nucCenterData, orbTypeData, orbExpData, 
						 molOrbiData,
						 rho_cut
						);
				

				//std::cout << "X1["<<IJK<<"]:" << x[IJK] << " Y1: " << y[IJK] << " Z1: "  << z[IJK] << " X2: " << X1 << " Y2: " << Y1 << " Z2: "  << Z1<<  " KLD: " << KLD << std::endl;

				}
			} else {
				KLD = 0.0;
			}
			
			
			long double valuePI = partial_int[III];
			long double w_beckeI = w_becke[IJK];
			partial_int[III] = valuePI + w_beckeI*KLD;
			
			//std::cout << "partial_int.at(" << III << ")= " << valuePI << std::endl;
			//std::cout << "+ w_becke[" << IJK <<"]:" << w_beckeI << "* KLD: " << KLD << std::endl;
			
			//std::cout << "partial_int.at(" << III << ")=" << partial_int.at(III) << std::endl;
		}
		
		
		double porcentaje = IJK;
		double count = int((porcentaje/npoints)*100);
		if (count>=progress){
			std::cout << "Processed:[" << IJK << "] of [" << npoints << "]:";
			std::cout << (porcentaje/npoints)*100 << " %" << std::endl;
			progress=progress+1;
		}
		

	}
	//write the results
	ofstream resultsFileOS;
	resultsFileOS.open (resultsFile);
	for (int III =0;III<IN1;III++){
		if (cubeXData[0]==0 && cubeYData[0] == 0){
			X1 = cubeZ + III*cubeZData[3];
		}
		if (cubeXData[0]==0 && cubeZData[0] == 0){
			X1 = cubeY + (III)*cubeYData[2];
		}
		if (cubeYData[0]==0 && cubeZData[0] == 0){
			X1 = cubeX + (III)*cubeXData[1];
		}
		resultsFileOS << "   " << X1 << "        " <<  partial_int[III] << "\n";
		
		//std::cout << "partial_int.at(" << III << ")="<<partial_int.at(III)<<endl;
	}
	resultsFileOS.close ();
	
	
	//releases reserved memory
	delete [] x;
	delete [] y;
	delete [] z;
	delete [] w_becke;
	delete [] partial_int;
	delete [] totalInt;
	
	return calc;
}

/**
 * reads input parameters from specified file when the program is executed
 */
int wfn_input (std::string paramFile){
    cout << "Reading file: "<< paramFile;
    std::ifstream file(paramFile.c_str());
    std::string str;
    
    //read initial parameter file: wfn.dat
    int line = 1;
    string dataFile;
    while (std::getline(file, str))
    {
        if (line == 1) {
            numberWFN = str;
            
        }
        if (line == 2){
            dataFile = str;
            //eliminates ' because it loads ' file name '
            dataFile = dataFile.substr(1,dataFile.length()-2);
            break;
        }
        line++;
    }
    file.close();
    cout << " OK"<<endl;
    
    //read the dataFile specified in the parameters by example water_dimer_rhf.wfn
    std::size_t found = paramFile.find_last_of("/\\");
    string dataFilePath = paramFile.substr(0,found+1);
    
    dataFilePath.append(dataFile);
    cout << "Reading data file:  "<< dataFilePath;
    std::ifstream dataPath(dataFilePath.c_str());
    line =1;
    
    /*
     * regular expressions to split the input line from the file
     */
    regex reg ("\\s+");
    bool loadNewRow = false;
    while (std::getline(dataPath, str))
    {
        //tokenized for str
        sregex_token_iterator iter(str.begin(), str.end(), reg, -1);
        sregex_token_iterator end;
        vector<string> fieldsStr (iter, end);
        
        //loads the parameters Gaussian, number of rows, number of columns, nucleis
        if (line == 2){
            method = fieldsStr[0];
            rowsMol = std::stoi(fieldsStr[1]);
            columsPrim = std::stoi(fieldsStr[4]);
            nuclei = std::stoi(fieldsStr[6]);
            cout << endl << "Method: "<< method << " Matrix: " << rowsMol << "x" << columsPrim << " and " << nuclei << endl;
        }
        //loads data according nuclei from the third line below n lines=nuclei
        /*  column 1  is the name
            column 4 is X
            column 5 is Y
            comumn 6 is Z
            column 7 is Charge
         */
        if (line >= 3 && line < nuclei+3){
	    //column 1: determines the position of the symbol in the array of Atoms
	    int symbolNumber = getSymbolNumber (fieldsStr[1]);
            rowNucleiData.push_back(symbolNumber);
            //column 4 to double
            std::istringstream os5(fieldsStr[5]);
            double number5;
            os5 >> number5;
            rowNucleiData.push_back(number5);
            //column 5 to double
            std::istringstream os6(fieldsStr[6]);
            double number6;
            os6 >> number6;
            rowNucleiData.push_back(number6);
            //column 6 to double
            std::istringstream os7(fieldsStr[7]);
            double number7;
            os7 >> number7;
            rowNucleiData.push_back(number7);
            //column 7 to double
            std::istringstream os10(fieldsStr[10]);
            double number10;
            os10 >> number10;
            rowNucleiData.push_back(number10);
            //puts the row and clear for next
            nucleiData.push_back(rowNucleiData);
            rowNucleiData.clear();
        }
        
        /**
         * LOADS NUCLEAR CENTER ASSOCIATED WITH EACH PRIMITIVE
         */ 
        if (str.compare(0,18,"CENTRE ASSIGNMENTS")==0){
            for (string a : fieldsStr)
            {
                if (a.compare (0,6,"CENTRE")!=0 && a.compare(0,11,"ASSIGNMENTS")!=0){
                    nucCenterData.push_back(stoi(a));
                }
            }
        }
        /**
         * LOADS TYPE OF ORBITAL FOR EACH PRIMITIVE
         */
        if (str.compare(0,16,"TYPE ASSIGNMENTS")==0){
            for (string a : fieldsStr)
            {
                if (a.compare (0,4,"TYPE")!=0 && a.compare(0,11,"ASSIGNMENTS")!=0){
                        orbTypeData.push_back(stoi(a));
                }
            }
        }
        /**
         * LOADS ORBITAL EXPONENTS FOR EACH PRIMITIVE
         */
        if (str.compare(0,9,"EXPONENTS")==0){
            for (string a : fieldsStr)
            {
                if (a.compare (0,9,"EXPONENTS")!=0){
                    //replace 0.000D+02 by 0.000e+02 to parsing according to sci notation
                    std::size_t pos = a.find("D");
                    if (pos != std::string::npos){
                        a.replace(pos, 1, "e");
                        std::istringstream os(a);
                        double number;
                        os >> number;
                        orbExpData.push_back(number);
                    }
                    
                }
            }
        }
        /**
         * LOADS MOLECULAR ORBITAL COEFFICIENTS
         */
        if (str.compare(0,2,"MO")==0 && loadNewRow == false){
            //the line begins with MO
            loadNewRow = true;
            //loads the occupational number of the line of MO (eight field)
            string a = fieldsStr[7];
            std::istringstream os(a);
            double number;
            os >> number;
            ocupNum.push_back(number);
            
        }else if (loadNewRow && str.compare(0,2,"MO")!=0){
            //loads each field of the line in the row after the line MO
            for (string a : fieldsStr)
            {
                //adds a new column in the row except the last data line
                if ((a.compare (0,3,"END")!=0 && (a.compare (0,4,"DATA"))!=0)){
                    //replace 0.000D+02 by 0.000e+02 to parsing according to sci notation
                    std::size_t pos = a.find("D");
                    if (pos != std::string::npos){
                        a.replace(pos, 1, "e");
                        std::istringstream os(a);
                        double number;
                        os >> number;
                        rowMolOrbiData.push_back(number);
                    }
                    
                }
            }
        }else if ((str.compare(0,2,"MO")==0 && loadNewRow)){
            //add the row to the molData when reaches new MO
            molOrbiData.push_back(rowMolOrbiData);
            //clear the row for new data
            rowMolOrbiData.clear();
            loadNewRow = true;
            //loads the occupational number of the line of MO (eight field)
            string a = fieldsStr[7];
            std::istringstream os(a);
            double number;
            os >> number;
            ocupNum.push_back(number);
        }
        if ((str.compare(0,3,"END")==0)){
            //add the last row to the molData when reaches END DATA
            molOrbiData.push_back(rowMolOrbiData);
            //clear the row for new data
            rowMolOrbiData.clear();
            loadNewRow = false;
        }
        line++;
    }
    dataPath.close();
    
    // validates the ocuppational number and number of electrons
    numberElectrons=0;
    for (double& fieldON : ocupNum){
        fieldON = fieldON * 0.5;
        if (fieldON>1) {
            fieldON=1;
        }else if (fieldON<0){
            fieldON=0;
        }
        numberElectrons=numberElectrons+fieldON;
        double elementAt = fieldON;
        if (fieldON==1){
            elementAt=0;
        }
        numberA.push_back(elementAt);
    }
        
    
    /**
     FOR DEBUG
     
    
    //prints nuclei data
    cout << endl << "Matrix of nuclei data is: " << endl;
    for (const auto& s1 : nucleiData) {
        cout << "BEGIN ROW"<< endl;
        for (const auto& s2 : s1)
            cout << s2 << " ";
        
        cout << "END ROW"<< endl;
    }
    //prints nuclear center data  of primitives
    cout << "Vector for nuclear  center of primitives is: " << endl;
    for (const auto& s1 : nucCenterData) {
        cout << s1 << " ";
    }
    
    //prints type of orbital data of primitives
    cout << endl << "Vector for type of orbital data for primitives is: " << endl;
    for (const auto& s1 : orbTypeData) {
        cout << s1 << " ";
    }
    
    //prints data of exponents for each primitive
    cout << endl << "Vector of orbital exponents is: " << endl;
     for (const auto& s1 : orbExpData) {
                 cout << s1 << " ";
     }
    //prints data of occupational number
    cout << endl << "Vector of occupational numbers is: " << endl;
    for (const auto& s1 : ocupNum) {
        cout << s1 << " ";
    }
    cout << "Number of electrons: "<< numberElectrons << endl;
    //prints data of A vector
    cout << endl << "Vector of A correction: " << endl;
    for (const double& s1 : numberA) {
        cout << s1 << " ";
    }
    //prints data of coefficients
    cout << endl << "Matrix orbital coefficients is: " << endl;
    for (const auto& s1 : molOrbiData) {
        cout << "BEGIN ROW"<< endl;
        for (const auto& s2 : s1)
            cout << s2 << " ";
        
        cout << "END ROW"<< endl;
    }
    */
	    
    cout << "  " <<  line << " lines.  OK"<<endl;
    
    // Executes Levenberg Marquardt Optimization for chiˆ2
    int result = optLevMarOpt ();
    
    /**
     * Reads cube.dat file
     */
    resultsFile = cubeInput (paramFile);
    /**
     * Becke Calculations
     */
    bool becke = becke88();
    
    return result;
}
/**
 validates if the file exists
 */
bool fileExists (const string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
}

/*
 * function to validate input parameters
 */
int paramsValidation (const char * app[], const char  * kldV[], const char * paramFile[]){
    int result=1;
    //checking correct versuion of kld
    if ((strcmp(kldV[0],"1D")==0)||(strcmp(kldV[0],"2D")==0)||(strcmp(kldV[0],"1D")==0)){
        result = 0;
    }else {
        cout << "Incorrect KLD version" << endl;
        result = 1;
    }
    //checking if the file exists
    if (fileExists (paramFile[0])){
        result =0;
    }else{
        cout << "Cant't open the specified file" << endl;
        result = 1;
    }
    
    return result;
}
/* Auxiliary routine: printing a matrix */
void print_matrix( string desc, int m, int n, double* a, int lda ) {
    int i, j;
    cout << desc ;
    for( i = 0; i < m; i++ ) {
        for( j = 0; j < n; j++ ) cout << a[i+j*lda];
        cout << endl;;
    }
}

/* Auxiliary routine: printing a vector of integers */
void print_int_vector( string desc, int n, int* a ) {
    int j;
    cout << desc;
    for( j = 0; j < n; j++ ) cout << a[j];
    cout << endl;;
}


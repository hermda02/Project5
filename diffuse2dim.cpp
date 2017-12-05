// Just beginning the final project for Computational Physics
// Diffusion Equation in 2 dimensions

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <armadillo>
#include <cstring>
#include <string>
#include <math.h>
#include <vector>
#include <fstream>
#include <time.h>
#include <ctime>

using namespace std;
using namespace arma;

ofstream ofile;

//  ----------------------------------------------//
//												  //
//	Two-dimensional diffusion equation solver:    //
//												  //
//	Utilizing armadillo library for automated     //
//	vector/matrix handling and memory allocation. //
//												  //
//--------------Matrix conditions-----------------//
//												  //
//		0	0	0	0	0	0	0	0	0		  //
//		0	0	0	0	0	0	0	0	0		  //
//		0	0	0	0	0	0	0	0	0		  //
//		0	0	0	0	0	0	0	0	0		  //
//		0	0	0	0	0	0	0	0	0		  //
//		0	0	0	0	0	0	0	0	0		  //
//		0	0	0	0	0	0	0	0	0		  //
//		0	0	0	0	0	0	0	0	0		  //
//		1	1	1	1	1	1	1	1	1		  //
//												  //
//  ----------------------------------------------//

//  Initializing the functions we will use:

//  Explicit Method
void ExplicitSolver(int,double,double,double,double);
//  Implicit Method
//void ImplicitSolver(int,double,double,double,double);
//  Tridiagonal Matrix Solver
void JacobiSolver(int,double,double,mat&,mat&,mat&,double);

//  ----------------------------Diffusion Equation Solver----------------------------

int main(int argc, char*argv[]){

	// Initializing the input arguments
	string filename;
	double tsize,latsize;
	int tsteps,style;

	if (argc <=4){
		cout << "Bad Usage:" << argv[0] << "read output file." << endl;
		exit(1);
	}

	// Saving input arguments
	if (argc > 1) {
		filename = argv[1];
		tsize = atof(argv[2]);
		latsize = atof(argv[3]);
		tsteps = atoi(argv[4]);
		style = atoi(argv[5]);
	}

	int n = 1/latsize;

	// To ensure explicit method stability
	if (tsize > (0.25*latsize*latsize)){
		tsize = 0.25*latsize*latsize;
	}

	//  --------------------------Initial Conditions---------------------------------
	double alpha = tsize/(latsize*latsize);
	double diag = 0;
	double offdiag = 0;
	double tolerance = 1.0e-8;

	// cout the inputs and their purposes for user verification
	cout << "File Name: " << filename << endl;
	cout << "Time step size: " << tsize << endl;
	cout << "lattice Step Size: " << latsize << endl;
	cout << "Time steps: " << tsteps << endl;

	// ------------------------Diffusion Equation Solvers----------------------------

	if (style == 1){
		cout << "Explicit Method" << endl;
		string fileout = filename;
		fileout += "_exp.txt";
		ofile.open(fileout);
		ofile << setiosflags(ios::showpoint | ios::uppercase);

		ExplicitSolver(n,alpha,tsize,tsteps,tolerance);
		ofile.close();
	}

	// if (style == 2){
	// 	cout << "Implicit Method" << endl;
	// 	string fileout = filename;
	// 	fileout += "_imp.txt";
	// 	ofile.open(fileout);
	// 	ofile << setiosflags(ios::showpoint | ios::uppercase);

	// 	ImplicitSolver(n,alpha,tsize,tsteps,tolerance);
	// 	ofile.close();
	//}

	if (style == 3){
		cout << "Doing Both Methods" << endl;

		// Explicit
		string fileout = filename;
		fileout += "_exp.txt";
		ofile.open(fileout);
		ofile << setiosflags(ios::showpoint | ios::uppercase);

		ExplicitSolver(n,alpha,tsize,tsteps,tolerance);
		ofile.close();

		// // Implicit
		// string fileout1 = filename;
		// fileout1 += "_imp.txt";
		// ofile.open(fileout1);
		// ofile << setiosflags(ios::showpoint | ios::uppercase);

		// ImplicitSolver(n,alpha,tsize,tsteps,tolerance);
		// ofile.close();
	}
}

void ExplicitSolver(int n, double alpha,double tsize, double tsteps,double tolerance){
	mat A = zeros<mat>(n,n);
	mat Aold = zeros<mat>(n,n);
	mat q = zeros<mat>(n,n);

	// Source term matrix -- heat comes from the bottom of the box
	for (int i=0;i<n;i++){
		q(n-1,i) = 1.0;
	}

	// cout << q << endl;

	for (int t =1;t<tsteps;t++){
		//  Initialize the "old" matrix. Make a guess at the values
		for (int i=0;i<n;i++){
			for (int j=0;j<n;j++){
				Aold(i,j) = 0.0;
			}
		}

		//  Boundary Contions, the lower boundary is equal to one.
		for (int i=0;i<n;i++){
			A(0,i) = 0.0;
			A(i,0) = 0.0;
			A(n-1,i) = 1.0;
			A(i,n-1) = 1.0;
		}
		JacobiSolver(n,tsize,alpha,A,Aold,q,tolerance);
	}
	for (int i=0;i<n;i++){
		for (int j=0;j<n;j++){
			ofile << setw(8) << setprecision(4) << A(i,j);
		}
		ofile << endl;
	}
	// cout << A << endl;
}

void JacobiSolver(int n,double dt,double alpha, mat &A, mat &Aold, mat &q, double abstol){
	int MaxIterations = 10000;
	
	//  Starting the iterative solver
	for (int k=0;k<MaxIterations;k++){
		for (int i=1;i<n-1;i++){
			for (int j=1;j<n-1;j++){
				A(i,j) = dt*q(i,j) + Aold(i,j) + 
				alpha*(Aold(i+1,j) + Aold(i,j+1) - 4.0*Aold(i,j) +
					   Aold(i-1,j) + Aold(i,j-1)); 
			}
		}
		double sum = 0.0;
		for (int i=0;i<n;i++){
			for (int j=0;j<n;j++){
				sum += (Aold(i,j)-A(i,j))*(Aold(i,j)-A(i,j));
				Aold(i,j) = A(i,j);
			}
		}
		if ( k == MaxIterations-1){
			if (sqrt(sum) >= abstol ){
			cout << "Matrix failed to converge!" << endl;
			}
		}
	}
}

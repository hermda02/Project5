// Just beginning the final project for Computational Physics
// Diffusion Equation

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
//	One-dimensional diffusion equation solver:    //
//												  //
//	Utilizing armadillo library for automated     //
//	vector/matrix handling and memory allocation. //
//												  //
//  ----------------------------------------------//


// Initialize the forward Euler matrix
void ForwardEuler(int,int,double,double);
// Initialize the backward Euler matrix
void BackwardEuler(int,int,double,double);
// Initialize the time centered matrix
void CrankNicol(int,int,double,double);

//Tridiagonal Solver
void TriSolver(int,double,double,vec,vec&);

// ----------------------------Diffusion Equation Solver-----------------------------

int main(int argc, char* argv[]){

	// Initialize the arguments defined by the inputs
	string filename;
	double tsize;
	double xsize;
	int tsteps,style;

	if (argc <= 4) {
		cout << "Bad Usage:" << argv[0] << "read output file ." << endl;
		exit(1);
	}

	// Saving arguments
	if (argc > 1) {
		filename=argv[1];
		tsize=atof(argv[2]);
		xsize=atof(argv[3]);
		tsteps=atoi(argv[4]);
		style=atoi(argv[5]);

	}

	// cout the inputs and their purposes for user verification
	cout << "File Name: " << filename << endl;
	cout << "Time step size: " << tsize << endl;
	cout << "X step size: " << xsize << endl;
	cout << "Time steps: " << tsteps << endl;
 
	// --------------Inputs Finished-----------------------
	// ------------Initial Conditions----------------------

	int n = 1/xsize;

	double alpha = tsize/(xsize*xsize);
	double diag =0;
	double offdia =0;

	// Diffusion Equation Solvers

	if (style == 1){
		cout << "Forward Euler Method" << endl;
		string fileout = filename;
		fileout += "_forward.txt";
		ofile.open(fileout);
		ofile << setiosflags(ios::showpoint | ios::uppercase);

		ForwardEuler(n,tsteps,xsize,alpha);
		ofile.close();
	}

	if (style == 2){
		cout << "Backward Euler Method" << endl;
		string fileout = filename;
		fileout += "_backward.txt";
		ofile.open(fileout);
		ofile << setiosflags(ios::showpoint | ios::uppercase);

		BackwardEuler(n,tsteps,xsize,alpha);
		ofile.close();
	}

	if (style == 3){
		cout << "Crank-Nicolson Method" << endl;
		string fileout = filename;
		fileout += "_crank.txt";
		ofile.open(fileout);
		ofile << setiosflags(ios::showpoint | ios::uppercase);

		CrankNicol(n,tsteps,xsize,alpha);
		ofile.close();
	}
	if (style == 4){
		cout << "Doing all three methods!" << endl;
		string fileout1 = filename;
		fileout1 += "_forward.txt";
		ofile.open(fileout1);
		ofile << setiosflags(ios::showpoint | ios::uppercase);

		ForwardEuler(n,tsteps,xsize,alpha);
		ofile.close();

		string fileout2 = filename;
		fileout2 += "_backward.txt";
		ofile.open(fileout2);
		ofile << setiosflags(ios::showpoint | ios::uppercase);

		BackwardEuler(n,tsteps,xsize,alpha);
		ofile.close();

		string fileout3 = filename;
		fileout3 += "_crank.txt";
		ofile.open(fileout3);
		ofile << setiosflags(ios::showpoint | ios::uppercase);

		CrankNicol(n,tsteps,xsize,alpha);
		ofile.close();
	}
}

void TriSolver(int n,double diag, double offdia, vec r, vec &u){
	// Initialize tridiagonal vectors (off diagonals are equal)
	// To be explicit, r is put in as the initial solution, and u is returned.
	// r is then redefined as the solution to the previous time step, and we repeat.
	// vector for new diagonal values
	vec c_star = zeros<vec>(n);
	vec u_star = zeros<vec>(n+1);

	// Forward Substitution
	c_star(0)=offdia/diag;
	u_star(0)=r(0)/diag;
	for (int i=1;i<n;i++){
		double m = 1.0/(diag-offdia*c_star(i-1));
		c_star(i)=offdia*m;
		u_star(i) = (r(i) - offdia*u_star(i-1))*m;
	}

	// Backward Substitution
	for (int i=n-1;i>=1;i--){
		u(i) = u_star(i)-c_star(i)*u(i+1);
	}
}

void ForwardEuler(int n, int tsteps, double xsize, double alpha){
	double diag, offdia;
	vec u = zeros<vec>(n+1);		// This is the vector u in Au = r
	vec unew = zeros<vec>(n+1);

	// Initial output file conditions
	for (int i=1;i<n;i++){
		double x = xsize*i;
		ofile << setw(6) << setprecision(2) << x;
	}
	ofile << endl;

	// Boundary Conditions
	unew(0) = u(0) = unew(n) = 0;
	u(n) = 1;
	diag = 1 - 2*alpha;
	offdia = alpha;

	// Time integration
	for (int t=1;t<=tsteps;t++){
		for (int i=1;i<n;i++){
			unew(i) = offdia*u(i-1) + diag*u(i) + offdia*u(i+1);
			u(i) = unew(i);
		}
		if (t % 5 ==0){
			for (int i=1;i<n;i++){
				ofile << setw(15) << setprecision(4) << u(i);
			}
			ofile << endl;
		}
		
	}
}

void BackwardEuler(int n, int tsteps, double xsize, double alpha){
	
	double diag, offdia;
	vec u = zeros<vec>(n+1);		// This is the vector u in Au = r
	vec r = zeros<vec>(n+1);		// The R.H.S. vector in Au = r

	// Initial vector conditions
	for (int i=1;i<n;i++){
		double x = xsize*i;
		ofile << setw(6) << setprecision(2) << x;
	}
	ofile << endl;

	// Boundary Conditions
	r(0) = u(0) = 0;
	r(n) = u(n) = 1;

	diag = 1 + 2*alpha;
	offdia = -1*alpha;
	// Time iteration here
	for (int t = 1; t <= tsteps; t++){
		// Using u as the previous time step solution
		TriSolver(n,diag,offdia,r,u);
		u(0) = 0;
		u(n) = 1;
		// replace previous time solution with new solution
		for (int i=0;i<=n;i++){
			r(i) = u(i);
		}
		if (t % 5 ==0){
			for (int i=1;i<n;i++){
				ofile << setw(15) << setprecision(4) << u(i);
			}
			ofile << endl;
		}
	}

}

void CrankNicol(int n, int tsteps, double xsize, double alpha){
	double diag, offdia;
	vec u = zeros<vec>(n+1);		// This is the vector u in Au = r
	vec unew = zeros<vec>(n+1);

	// Initial output file conditions
	for (int i=1;i<n;i++){
		double x = xsize*i;
		ofile << setw(6) << setprecision(2) << x;
	}
	ofile << endl;

	// Boundary Conditions
	unew(0) = u(0) = 0;
	unew(n) = u(n) = 1;
	diag = 1 - alpha;
	offdia = 0.5*alpha;

	for (int t=1;t <= tsteps; t++){
		for (int i=1;i<n;i++){
			unew(i) = offdia*u(i-1) + diag*u(i) + offdia*u(i+1);
		}
		for (int i=1;i<n;i++){
			u(i) = unew(i);
		}
		//TriSolver(n,diag,offdia,alpha,u);
		if (t % 5 ==0){
			for (int i=1;i<n;i++){
				ofile << setw(15) << setprecision(4) << u(i);
			}
			ofile << endl;
		}
	}
}
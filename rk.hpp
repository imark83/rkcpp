
#ifndef __RK_HPP__
#define __RK_HPP__
// performs an RK step. Returns -1 if rejected,
// otherwise returns FACTOR for next step size

// typedef void (*cardioFun)(double, double*, double*, double*);

#include "variable.hpp"

#define MAX(i,j) ( ((i) > (j)) ? (i) : (j) )
#define MIN(i,j) ( ((i) > (j)) ? (j) : (i) )

void computeStages (int nvar, 		// number of variables
			Variable rkStage[],					// RK stages
			Variable x[],								// integrated variable
			double t,										// integration variable
	 		double h, 									// step size
			double *pars);							// parameters
			// cardioFun f);

//
void rk(int nvar, 						// number of variables of dependent variable
	Variable x[], 							// dependent variable
	double t0, 									// initial time
	double tf, 									// end time
	double denseStep,						// dense step for output
	double *pars, 							// parameters
	double tol,									// parameters
	int event,									// variable to compute poincare sections. -1 none
	double eventVal);						// poincare section value

#endif

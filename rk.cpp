#include <iostream>
#include <fstream>
#include <cstdio>
#include <math.h>
#include "rk.hpp"
#include "rk_coefs.hpp"
#include "fun.hpp"

using namespace std;

// extern int nsteps;
// extern int nrejected;



void computeStages(int nvar, 		// number of variables
			double rkStage[],					// RK stages
			double x[],								// integrated variable
			double t,										// integration variable
	 		double h, 									// step size
			double *pars) {							// parameters
			// cardioFun f);

	// COMPUTE STAGES OF THE STEP(assume we have an autonomous system)
	// 2
	for(int j=0; j<nvar; ++j) rkStage[1*nvar+j] = x[j]
					+ h *(A21*rkStage[j]);
	f(t + C2*h,(rkStage + 1*nvar),(rkStage + 1*nvar), pars);
	// 3
	for(int j=0; j<nvar; ++j) rkStage[2*nvar+j] = x[j]
					+ h *(A31*rkStage[j] + A32*rkStage[1*nvar+j]);
	f(t + C3*h,(rkStage + 2*nvar),(rkStage + 2*nvar), pars);
	// 4
	for(int j=0; j<nvar; ++j) rkStage[3*nvar+j] = x[j]
					+ h *(A41*rkStage[j] + A42*rkStage[1*nvar+j]
					+ A43*rkStage[2*nvar+j]);
	f(t + C4*h,(rkStage + 3*nvar),(rkStage + 3*nvar), pars);
	// 5
	for(int j=0; j<nvar; ++j) rkStage[4*nvar+j] = x[j]
					+ h *(A51*rkStage[j] + A52*rkStage[1*nvar+j]
					+ A53*rkStage[2*nvar+j] + A54*rkStage[3*nvar+j]);
	f(t + C5*h,(rkStage + 4*nvar),(rkStage + 4*nvar), pars);
	// 6
	for(int j=0; j<nvar; ++j) rkStage[5*nvar+j] = x[j]
					+ h *(A61*rkStage[j] + A62*rkStage[1*nvar+j]
					+ A63*rkStage[2*nvar+j] + A64*rkStage[3*nvar+j]
					+ A65*rkStage[4*nvar+j]);
	f(t + C6*h,(rkStage + 5*nvar),(rkStage + 5*nvar), pars);
	// 7(A72 = 0)
	for(int j=0; j<nvar; ++j) rkStage[6*nvar+j] = x[j]
					+ h *(A71*rkStage[j] + A73*rkStage[2*nvar+j]
					+ A74*rkStage[3*nvar+j] + A75*rkStage[4*nvar+j]
					+ A76*rkStage[5*nvar+j]);
	f(t + C7*h,(rkStage + 6*nvar),(rkStage + 6*nvar), pars);

	return;
}



double estimateError(int nvar, double rkStage[], double h) {
	// ESTIMATE ERROR(E2 = 0)
	double error = 0.0;
	for(int j=0; j<nvar; ++j) error = error + fabs(E1*rkStage[j] + E3*rkStage[2*nvar+j] + E4*rkStage[3*nvar+j] + E5*rkStage[4*nvar+j] + E6*rkStage[5*nvar+j] + E7*rkStage[6*nvar+j]);
	return error * h;
}



int spikePoincare(int nvar, double t, double step, double x[],
				double xNext[], double rkStage[], double tol, double pars[],
				double &eventT, double eventX[]) {


	double fx_L[nvar], fx_M[nvar], fx_R[nvar];
	double th_L = 0.0;
	double th_R = 1.0;
	double th_M = 0.5;
	for(int j=0; j<nvar; ++j){
		fx_L[j] = rkStage[j];
	}
	f(t, fx_R, xNext, pars);

	// double eventX[nvar];
	if(fx_L[0] > 0 && fx_R[0] < 0) {
		// THERE IS A MAXIMUM(SPIKE)
		for(int j=0; j<nvar; ++j)
			fx_M[j] = denseEval(nvar, rkStage, x, step, j, th_M);
		f(t, fx_M, fx_M, pars);
		do {
			if(fx_M[0] > 0.0) {
				th_L = th_M;
				th_M = (th_R + th_L) / 2.0;
				for(int j=0; j<nvar; ++j) {
					fx_L[j] = fx_M[j];
					fx_M[j] = denseEval(nvar, rkStage, x, step, j, th_M);
				}
				f(t, fx_M, fx_M, pars);
			} else {
				th_R = th_M;
				th_M = (th_R + th_L) / 2.0;
				for(int j=0; j<nvar; ++j) {
					fx_R[j] = fx_M[j];
					fx_M[j] = denseEval(nvar, rkStage, x, step, j, th_M);
				}
				f(t, fx_M, fx_M, pars);
			}
		} while(fabs(fx_M[0]) > 10.0*tol);
		eventT = t+th_M * step;
		for(int j=0; j<nvar; ++j) {
			eventX[j] = denseEval(nvar, rkStage, x, step, j, th_M);
		}
		return 1;
	}
	return 0;
}


//
void rk(int nvar, 			// number of variables of dependent variable
			double x[], 				   	// dependent variable
			double t0, 						  // initial time
			double tf, 							// end time
			double denseStep,				// dense step for output
			double *pars, 					// parameters
			double tol,							// parameters
			int event,							// 0 -> none. 1-> period. 2->nspikes and DC
			Task &task) {


	double step = 1.0e-6;			// INTEGRATION STEP SIZE
	double rkStage[7*nvar];		// RK STAGES
	double t = t0;							// integration time
	double denseT = t0;							// time for dense output



	// VARIABLES FOR SPIKE NUMBER COMPUTATION
	double xNext[nvar];

	const int MAX_SPIKES = 64;
	char nSpikes = 0;
	double refSpike[nvar];
	double eventT[MAX_SPIKES], eventX[nvar];


	// INITIALIZE FSAL STAGE(store in stage 0)
	for(int j=0; j<nvar; ++j) rkStage[j] = x[j];

	f(t, rkStage, rkStage, pars);

	double fac;									// step size correction factor
	char endOfIntegration = 0;	// end of integration flag


   std::cout << denseT;
   for(int j=0; j<3; ++j)
       std::cout << "  " << x[j];
   std::cout << std::endl;

	// MAIN LOOP
	while(!endOfIntegration) {
		// ++nsteps;
		computeStages(nvar, rkStage, x, t, step, pars);
		double error = estimateError(nvar, rkStage, step);


		// compute norm of x
		double normX = fabs(x[0]);
		for(int j=1; j<nvar; ++j) normX += fabs(x[j]);

		// CHECK IF WE REJECT STEP
		if(error > tol || error > tol*normX) {
			// ++nrejected;
			step = 0.2 * step;
			continue;
		}

		// UP TO HERE ACCEPTED STEP
		if(tf-t < step) {
			endOfIntegration = 1;
			step = tf-t;
      computeStages(nvar, rkStage, x, t, step, pars);
		}


		// USE THE 5TH ORDER RK TO GO AHEAD(B2 = B7 = 0)
		for(int j=0; j<nvar; ++j) xNext[j] = x[j] + step *(B1*rkStage[j] +
			 		B3*rkStage[2*nvar+j] + B4*rkStage[3*nvar+j] + B5*rkStage[4*nvar+j] +
					B6*rkStage[5*nvar+j]);

		// POINCARE SECTION
    if(event==1) {
      if(spikePoincare(nvar, t, step, x, xNext, rkStage, tol, pars,
							eventT[nSpikes], eventX)) {
				// SPIKE FOUND
				if(nSpikes == 0) { 	// FIRST SPIKE
					++nSpikes;
					for(int j=0; j<nvar; ++j) {
						refSpike[j] = eventX[j];
					}
				} else {
					normX = fabs(refSpike[0]-eventX[0]);
					for(int j=1; j<nvar; ++j){
						if (fabs(refSpike[j] - eventX[j]) > normX) {
							normX = fabs(refSpike[j] - eventX[j]);
						}
					}
					if(normX < 1.0e-4) {

						// LOOP COMPLETE
						task.result.period = eventT[1] - eventT[0];
						return;
					}
				}
			}
    }
    if(event==2) { // COMPUTE SPIKES UP TO A MAXIMUM OF 64
      if(spikePoincare(nvar, t, step, x, xNext, rkStage, tol, pars,
							eventT[nSpikes], eventX)) {
				// SPIKE FOUND
				if(nSpikes == 0) { 	// FIRST SPIKE
					for(int j=0; j<nvar; ++j) {
						refSpike[j] = eventX[j];
					}
				} else {
					normX = fabs(refSpike[0]-eventX[0]);
					for(int j=1; j<nvar; ++j){
						if (fabs(refSpike[j] - eventX[j]) > normX) {
							normX = fabs(refSpike[j] - eventX[j]);
						}
					}
					if(normX < 1.0e-4 || nSpikes == MAX_SPIKES) {
						// LOOP COMPLETE
						task.result.sn |= (1 << nSpikes);
						double relax = eventT[1] - eventT[0];
						for(int i=2; i<nSpikes; ++i) {
							if(eventT[i] - eventT[i-1] > relax)
								relax = eventT[i] - eventT[i-1];
            }
            task.result.period = eventT[nSpikes] - eventT[0];
						task.result.dutyCycle = task.result.period - relax;
						return;
					}
				}
				++nSpikes;
				if(nSpikes == MAX_SPIKES) {
					task.result.sn |= 0x8000000000000000;
					task.result.dutyCycle = task.result.period;
				}
			}
    }

    // std::cout << x[1] << std::endl;

		// DENSE OUTPUT
		while(denseT + denseStep - t - step < 1.0e-15) {
			denseT += denseStep;
			double th = (denseT - t) / step;
           std::cout << denseT;
           for(int j=0; j<3; ++j)
               std::cout << "  " << denseEval(nvar, rkStage, x, step, j,th);
           std::cout << std::endl;
		}


		// USE THE 5TH ORDER RK TO GO AHEAD(B2 = B7 = 0)
		for(int j=0; j<nvar; ++j) x[j] = xNext[j];
		t += step;



		// ESTIMATE FACTOR FOR NEXT STEP SIZE
		fac = 0.8*pow((tol/error), 0.2);	// multiplier for next step size
		fac = MAX(fac, 0.1);			// no smaller than 10%
		fac = MIN(fac, 1.5);			// no bigger than 150%


		step = fac * step;

		// UPDATE FSAL
		for(int j=0; j<nvar; ++j) rkStage[j] = rkStage[6*nvar+j];

	}
	return;
}

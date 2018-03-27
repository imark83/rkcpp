#include <iostream>




// SECOND MEMBER FOR HENON-HEILES HAMILTONIAN SYSTEM
void f(double t, double *rop, double *op, double *p) {
	int nvar = 4;
	double x[nvar];
	for(int i=0; i<nvar; ++i) x[i] = op[i];

	// rop[0] = x[1];
	// rop[1] = -x[0];

	rop[0] = x[2];
	rop[1] = x[3];
	rop[2] = -1.0 * (x[0] + 2.0 * x[0] * x[1]);
	rop[3] = x[1] * x[1] - x[1] - x[0] * x[0];


  // rop[0] = 10.0 * (x[1] - x[0]);
  // rop[1] = x[0] * (28.0 - x[2]) - x[1];
  // rop[2] = x[0] * x[1] - 8.0/3.0 * x[2];


	return;
}

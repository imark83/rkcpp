#include <iostream>
#include <cmath>



// SECOND MEMBER FOR RETARDED COPULED NEURONS

void f(double t, double *rop, double *x, double *p) {
  int nvar=3;
  double aux[nvar];


  // p[0] -> b
  // p[1] -> I
  const double eps = 0.01;
  const double x0 = -1.6;

  // double GV;

  // Neurona 0

  aux[0] = x[1] + x[0]*x[0]*(p[0]-x[0]) - x[2] + p[1];
  aux[1] = 1.0 - 5.0*x[0]*x[0] - x[1];
  aux[2] = eps*(4.0*(x[0]-x0) - x[2]);


  rop[0] = aux[0];
  rop[1] = aux[1];
  rop[2] = aux[2];

	return;
}

#include <iostream>
#include <cmath>
#include "buffer.hpp"



// SECOND MEMBER FOR RETARDED COPULED NEURONS

void f(double t, double *rop, double *x, double *p) {
  int nvar=3;
  double aux[nvar];

  // p[0] -> vthKS
  // p[1] -> Iext

  double ninf, minf, cinf, taum, tauc, ICa, IK, IKS, IL, GV;

  // Neurona 0
  ninf = 1/(1+exp(-2*(0.055)*(x[0]-(-1.2))));
  minf = 1/(1+exp(-2*(0.1)*(x[0]-(2.0))));
  cinf = 1/(1+exp(-2*(0.4)*(x[0]-(p[0]))));

  taum = cosh ((0.1) * (x[0] - (2.0))/2.0);
  tauc = cosh ((0.4) * (x[0] - (p[0]))/2.0);

  ICa = (4.4) * ninf * (x[0]-(120.0));
  IK = (8.0) * x[1] * (x[0]-(-80.0));
  IKS = (0.15) * x[2] * (x[0]-(-80.0));
  IL = (2.0) * (x[0]-(-60.0));

  aux[0] = (-(ICa + IK + IL + IKS) + (p[1]))/(1.2);
  aux[1] = ((4.9) * (minf - x[1])) * taum;
  aux[2] = ((0.005) * (cinf - x[2])) * tauc;


  rop[0] = aux[0];
  rop[0] = aux[1];
  rop[2] = aux[2];

	return;
}

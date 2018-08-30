#include <iostream>
#include <cmath>



double sigmoid(double x) {
  return 1.0/(1.0 + exp(80*(x)));
}

// SECOND MEMBER FOR RETARDED COPULED NEURONS
void f(double t, double *rop, double *x, double *p) {

  // p[0] -> b
  // p[1] -> Iext
  // p[2] -> g_syn

  int nvar=18;

  double aux[nvar];


  // p[0] -> b
  // p[1] -> I
  const double eps = 0.01;
  const double x0 = -1.6;

  // double GV;

  // Neurona 0
  for(int k=0; k<6; ++k) {
    aux[3*k+0] = x[3*k+1] + x[3*k+0]*x[3*k+0]*(p[0]-x[3*k+0])
              - x[3*k+2] + p[1];
    aux[3*k+1] = 1.0 - 5.0*x[3*k+0]*x[3*k+0] - x[3*k+1];
    aux[3*k+2] = eps*(4.0*(x[3*k+0]-x0) - x[3*k+2]);
  }


  // COUPLING

  const double EPOST = -2.0;
  const double ETH = -1.2;

  aux[ 0] += p[2]*(x[ 0]-EPOST) * (sigmoid(x[3]-ETH) + sigmoid(x[9]-ETH));
  aux[ 3] += p[2]*(x[ 3]-EPOST) * (0.5*sigmoid(x[0]-ETH)
                  + 0.5*sigmoid(x[6]-ETH) + sigmoid(x[12]-ETH));
  aux[ 6] += p[2]*(x[ 6]-EPOST) * (sigmoid(x[3]-ETH) + sigmoid(x[15]-ETH));
  aux[ 9] += p[2]*(x[ 9]-EPOST) * (sigmoid(x[0]-ETH) + sigmoid(x[12]-ETH));
  aux[12] += p[2]*(x[12]-EPOST) * (0.5*sigmoid(x[9]-ETH) + sigmoid(x[3]-ETH)
                  + 0.5*sigmoid(x[15]-ETH));
  aux[15] += p[2]*(x[15]-EPOST) * (sigmoid(x[12]-ETH) + sigmoid(x[6]-ETH));


  for(int k=0; k<nvar; ++k)
    rop[k] = aux[k];
}

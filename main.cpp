#include <iostream>
#include <cmath>
#include "rk.hpp"
#include "buffer.hpp"



int nsteps = 0;
int nrejected = 0;
Buffer retard0(50000, 0.0);
Buffer retard1(50000, 0.0);
Buffer retard2(50000, 0.0);

int main(int argc, char const *argv[]) {
  std::cout.precision(15);
  std::cout << std::scientific;

  int nvar = 12;
  double x[nvar];

  x[0] = -1.0;
  x[1] = 0.0;
  x[2] = 0.0;
  x[3] = -1.5;
  x[4] = 0.0;
  x[5] = 0.0;
  x[6] = -2.0;
  x[7] = 0.0;
  x[8] = 0.0;
  x[9] = 0.0;
  x[10] = 0.0;
  x[11] = 0.0;


  double pars[2] = {30.0, -25.111111111111};


  double poincareThresHold = -30.0;


  rk(nvar, x, 0, 20000, 0.05, pars, 1.0e-8, -1, poincareThresHold);

  std::cerr << "nsteps = " << nsteps << std::endl;
  std::cerr << "nrejected = " << nrejected << std::endl;


  return 0;
}

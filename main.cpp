#include <iostream>
#include "rk.hpp"
#include <cmath>



int nsteps = 0;
int nrejected = 0;

int main(int argc, char const *argv[]) {
  std::cout.precision(15);
  std::cout << std::scientific;

  int nvar = 4;
  double x[nvar];

  x[0] = 0.1; x[1] = 0.0;
  x[2] = 0.0; x[3] = 0.05;

  rk(nvar, x, 0, 256, 256, NULL, 1.0e-8, -1, 0.0);

  std::cerr << "nsteps = " << nsteps << std::endl;
  std::cerr << "nrejected = " << nrejected << std::endl;


  return 0;
}

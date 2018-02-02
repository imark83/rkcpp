#include "variable.hpp"
#include "rk.hpp"
#include <cmath>



int nsteps = 0;
int nrejected = 0;

int main(int argc, char const *argv[]) {
  std::cout.precision(15);
  std::cout << std::scientific;

  int nvar = 2;
  Variable x[nvar];

  x[0] = 1.0; x[1] = 0;

  rk(nvar, x, 0, 10*M_PI, 10*M_PI, NULL, 1.0e-8);

  std::cerr << "nsteps = " << nsteps << std::endl;
  std::cerr << "nrejected = " << nrejected << std::endl;


  return 0;
}

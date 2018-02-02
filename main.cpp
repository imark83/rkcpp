#include "variable.hpp"
#include "rk.hpp"

int nsteps = 0;
int nrejected = 0;



int main(int argc, char const *argv[]) {
  std::cout.precision(15);

  int nvar = 4;
  Variable x[nvar];

  x[0] = 0.01; x[1] = 0; x[2] = 0; x[3] = 0.1;


  rk(nvar, x, 0, 500, 1, NULL, 1.0e-8);




  return 0;
}

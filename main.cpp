#include <iostream>
#include <fstream>
#include <cmath>
#include "rk.hpp"
#include "buffer.hpp"



int nsteps = 0;
int nrejected = 0;
Buffer retard0(50000, 0.0);
Buffer retard1(50000, 0.0);
Buffer retard2(50000, 0.0);
std::fstream fPoinc;

int main(int argc, char const *argv[]) {
  std::cout.precision(8);
  std::cout << std::scientific;

  int nvar = 12;
  double x[nvar];

  x[0] = -1.0;
  x[1] = 0.0;
  x[2] = 0.0;
  x[3] = -1.5;
  x[4] = 0.0;
  x[5] = 0.0;
  x[6] = -20.0;
  x[7] = 0.0;
  x[8] = 0.0;
  x[9] = 0.0;
  x[10] = 0.0;
  x[11] = 0.0;


  double pars[3] = {75, -26.77777777777, 0.005};


  double poincareThresHold = -30.0;
  fPoinc.open("spider.txt", std::fstream::app);

  rk(nvar, x, 0, 40000, 40000, pars, 1.0e-8, 1, poincareThresHold);
  fPoinc << std::endl;
  fPoinc.close();
  std::cerr << "nsteps = " << nsteps << std::endl;
  std::cerr << "nrejected = " << nrejected << std::endl;


  return 0;
}

#include <iostream>
#include <fstream>
#include <cmath>
#include "rk.hpp"



int nsteps = 0;
int nrejected = 0;

std::fstream fPoinc;

int main(int argc, char const *argv[]) {
  std::cout.precision(8);
  std::cout << std::scientific;

  int nvar = 12;
  double y[nvar];
  // Mierdas de Buffers
  Buffer retard[] = {Buffer(50000, 0.0), Buffer(50000, 0.0), Buffer(50000, 0.0)};
  Buffer auxBuffer, canonicalBuffer;

  y[0] = -1.0;
  y[1] = 0.0;
  y[2] = 0.0;
  y[3] = -1.5;
  y[4] = 0.0;
  y[5] = 0.0;
  y[6] = -20.0;
  y[7] = 0.0;
  y[8] = 0.0;
  y[9] = 0.0;
  y[10] = 0.0;
  y[11] = 0.0;


  // START WITH DECOUPLED NETWORK
  double pars[3] = {75, -26.77777777777, 0.0};
  double poincareThresHold = -30.0;
  rk(nvar, y, 0.0, 4000, 4000,
        pars, 1.0e-8, 0, poincareThresHold, retard);
  retard[0].resetTime();
  retard[1].resetTime();
  retard[2].resetTime();

  double P = rk(nvar, y, 0.0, 1000, 1000,
        pars, 1.0e-8, 1, poincareThresHold, retard);
  retard[0].resetTime();
  retard[1].resetTime();
  retard[2].resetTime();
  std::cout << "P = " << P << std::endl;
  canonicalBuffer = retard[0];


  double x[nvar];     // initial conditions for network
  double z[nvar];     // variables to rock & roll the neuron

  for(int i=0; i<3; ++i) {
    x[i] = z[i] = y[i];
  }
  x[9] = z[9] = y[9];
  rk(nvar, z, 0.0, (P*(1.0-0.45)), (P*(1.0-0.15)),
      pars, 1.0e-8, 0, poincareThresHold, retard);
  retard[0].resetTime();
  retard[1].resetTime();
  retard[2].resetTime();
  for(int i=0; i<3; ++i)
    x[3+i] = z[i];
  auxBuffer = retard[0];
  x[10] = z[9];


  // GRIND IT... AGAIN!
  for(int i=0; i<3; ++i) {
    z[i] = y[i];
  }
  retard[0] = canonicalBuffer;
  rk(nvar, z, 0.0, (P*(1.0-0.55)), (P*(1.0-0.85)),
      pars, 1.0e-8, 0, poincareThresHold, retard);
  retard[0].resetTime();
  retard[1].resetTime();
  retard[2].resetTime();

  for(int i=0; i<3; ++i)
    x[6+i] = z[i];
  x[11] = z[9];
  retard[2] = retard[0];
  retard[1] = auxBuffer;
  retard[0] = canonicalBuffer;


  pars[2] = 0.004;
  fPoinc.open("spider.txt", std::fstream::app);
  rk(nvar, x, 0.0, 5000, 0.05, pars,
        1.0e-8, 0, poincareThresHold, retard);
  fPoinc << std::endl;
  fPoinc.close();
  std::cerr << "nsteps = " << nsteps << std::endl;
  std::cerr << "nrejected = " << nrejected << std::endl;


  return 0;
}

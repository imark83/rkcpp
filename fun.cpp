#include <iostream>
#include <cmath>
#include "buffer.hpp"

extern Buffer retard0;
extern Buffer retard1;
extern Buffer retard2;


// SECOND MEMBER FOR RETARDED COPULED NEURONS

void f(double t, double *rop, double *x, double *p) {

  double aux[12];

  // p[0] -> DELAY
  // p[1] -> vthKS
  // p[2] -> g_syn

  double ninf, minf, cinf, taum, tauc, ICa, IK, IKS, IL, GV;

  // Neurona 0
  ninf = 1/(1+exp(-2*(0.055)*(x[0]-(-1.2))));
  minf = 1/(1+exp(-2*(0.1)*(x[0]-(2.0))));
  cinf = 1/(1+exp(-2*(0.4)*(x[0]-(p[1]))));

  taum = cosh ((0.1) * (x[0] - (2.0))/2.0);
  tauc = cosh ((0.4) * (x[0] - (p[1]))/2.0);

  ICa = (4.4) * ninf * (x[0]-(120.0));
  IK = (8.0) * x[1] * (x[0]-(-80.0));
  IKS = (0.15) * x[2] * (x[0]-(-80.0));
  IL = (2.0) * (x[0]-(-60.0));

  aux[0] = (-(ICa + IK + IL + IKS) + (35.5))/(1.2);
  aux[1] = ((4.9) * (minf - x[1])) * taum;
  aux[2] = ((0.005) * (cinf - x[2])) * tauc;

  // Neurona 1
  ninf = 1/(1+exp(-2*(0.055)*(x[3]-(-1.2))));
  minf = 1/(1+exp(-2*(0.1)*(x[3]-(2.0))));
  cinf = 1/(1+exp(-2*(0.4)*(x[3]-(p[1]))));

  taum = cosh ((0.1) * (x[3] - (2.0))/2.0);
  tauc = cosh ((0.4) * (x[3] - (p[1]))/2.0);

  ICa = (4.4) * ninf * (x[3]-(120.0));
  IK = (8.0) * x[4] * (x[3]-(-80.0));
  IKS = (0.15) * x[5] * (x[3]-(-80.0));
  IL = (2.0) * (x[3]-(-60.0));

  aux[3] = (-(ICa + IK + IL + IKS) + (35.5))/(1.2);
  aux[4] = ((4.9) * (minf - x[4])) * taum;
  aux[5] = ((0.005) * (cinf - x[5])) * tauc;

  // Neurona 2
  ninf = 1/(1+exp(-2*(0.055)*(x[6]-(-1.2))));
  minf = 1/(1+exp(-2*(0.1)*(x[6]-(2.0))));
  cinf = 1/(1+exp(-2*(0.4)*(x[6]-(p[1]))));

  taum = cosh ((0.1) * (x[6] - (2.0))/2.0);
  tauc = cosh ((0.4) * (x[6] - (p[1]))/2.0);

  ICa = (4.4) * ninf * (x[6]-(120.0));
  IK = (8.0) * x[7] * (x[6]-(-80.0));
  IKS = (0.15) * x[8] * (x[6]-(-80.0));
  IL = (2.0) * (x[6]-(-60.0));

  aux[6] = (-(ICa + IK + IL + IKS) + (35.5))/(1.2);
  aux[7] = ((4.9) * (minf - x[7])) * taum;
  aux[8] = ((0.005) * (cinf - x[8])) * tauc;

  // Synapsis 1->0
  aux[0] += -(p[2]) * x[10] * (x[0] - (-70.0)) / (1.2);
  // Synapsis 1->2
  aux[6] += -(p[2]) * x[10] * (x[6] - (-70.0)) / (1.2);
  // Sinapsis 0->1
  aux[3] += -(p[2]*0.5) * x[9] * (x[3] - (-70.0)) / (1.2);
  // Sinapsis 2->1
  aux[3] += -(p[2]*0.5) * x[11] * (x[3] - (-70.0)) / (1.2);

  // Synapsis pasada
  aux[0] += -(p[2]) * retard0(t-p[0]) * (x[0] - (-70.0)) / (1.2);
  aux[3] += -(p[2]) * retard1(t-p[0]) * (x[3] - (-70.0)) / (1.2);
  aux[6] += -(p[2]) * retard2(t-p[0]) * (x[6] - (-70.0)) / (1.2);

  // Actualizar variable sinaptica de N0
  GV = (0.002) / (1 + exp(-(0.22) * (x[0] - (2.0))));
  aux[9] = (5000.0) * GV * (1-x[9]) - (0.18)*x[9];

  // Actualizar variable sinaptica de N1
  GV = (0.002) / (1 + exp(-(0.22) * (x[3] - (2.0))));
  aux[10] = (5000.0) * GV * (1-x[10]) - (0.18)*x[10];
  // Actualizar variable sinaptica de N2
  GV = (0.002) / (1 + exp(-(0.22) * (x[6] - (2.0))));
  aux[11] = (5000.0) * GV * (1-x[11]) - (0.18)*x[11];

  rop[0] = aux[0];
  rop[1] = aux[1];
  rop[2] = aux[2];
  rop[3] = aux[3];
  rop[4] = aux[4];
  rop[5] = aux[5];
  rop[6] = aux[6];
  rop[7] = aux[7];
  rop[8] = aux[8];
  rop[9] = aux[9];
  rop[10] = aux[10];
  rop[11] = aux[11];

	return;
}

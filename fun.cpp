#include <iostream>
#include <cmath>
#include "buffer.hpp"



// SECOND MEMBER FOR RETARDED COPULED NEURONS
void f(double t, double *rop, double *x, double *p) {
  // p[0] -> vthKS
  // p[1] -> g_syn
  double aux[24];
  {
      double ninf, minf, cinf, taum, tauc, ICa, IK, IKS, IL;

      ninf = 1/(1+exp(-2*(0.055)*(x[3]-(-1.2))));
      minf = 1/(1+exp(-2*(0.1)*(x[3]-(2.0))));
      cinf = 1/(1+exp(-2*(0.4)*(x[3]-(p[0]))));

      taum = cosh ((0.1) * (x[3] - (2.0))/2.0);
      tauc = cosh ((0.4) * (x[3] - (p[0]))/2.0);

      ICa = (4.4) * ninf * (x[3]-(120.0));
      IK = (8.0) * x[4] * (x[3]-(-80.0));
      IKS = (0.15) * x[5] * (x[3]-(-80.0));
      IL = (2.0) * (x[3]-(-60.0));

      aux[3] = (-(ICa + IK + IL + IKS) + (36.0))/(1.2);
      aux[4] = ((4.9) * (minf - x[4])) * taum;
      aux[5] = ((0.005) * (cinf - x[5])) * tauc;
  }
  {
      double ninf, minf, cinf, taum, tauc, ICa, IK, IKS, IL;

      ninf = 1/(1+exp(-2*(0.055)*(x[15]-(-1.2))));
      minf = 1/(1+exp(-2*(0.1)*(x[15]-(2.0))));
      cinf = 1/(1+exp(-2*(0.4)*(x[15]-(p[0]))));

      taum = cosh ((0.1) * (x[15] - (2.0))/2.0);
      tauc = cosh ((0.4) * (x[15] - (p[0]))/2.0);

      ICa = (4.4) * ninf * (x[15]-(120.0));
      IK = (8.0) * x[16] * (x[15]-(-80.0));
      IKS = (0.15) * x[17] * (x[15]-(-80.0));
      IL = (2.0) * (x[15]-(-60.0));

      aux[15] = (-(ICa + IK + IL + IKS) + (36.0))/(1.2);
      aux[16] = ((4.9) * (minf - x[16])) * taum;
      aux[17] = ((0.005) * (cinf - x[17])) * tauc;
  }
  {
      double ninf, minf, cinf, taum, tauc, ICa, IK, IKS, IL;

      ninf = 1/(1+exp(-2*(0.055)*(x[0]-(-1.2))));
      minf = 1/(1+exp(-2*(0.1)*(x[0]-(2.0))));
      cinf = 1/(1+exp(-2*(0.4)*(x[0]-(p[0]))));

      taum = cosh ((0.1) * (x[0] - (2.0))/2.0);
      tauc = cosh ((0.4) * (x[0] - (p[0]))/2.0);

      ICa = (4.4) * ninf * (x[0]-(120.0));
      IK = (8.0) * x[1] * (x[0]-(-80.0));
      IKS = (0.15) * x[2] * (x[0]-(-80.0));
      IL = (2.0) * (x[0]-(-60.0));

      aux[0] = (-(ICa + IK + IL + IKS) + (36.0))/(1.2);
      aux[1] = ((4.9) * (minf - x[1])) * taum;
      aux[2] = ((0.005) * (cinf - x[2])) * tauc;
  }
  {
      double ninf, minf, cinf, taum, tauc, ICa, IK, IKS, IL;

      ninf = 1/(1+exp(-2*(0.055)*(x[9]-(-1.2))));
      minf = 1/(1+exp(-2*(0.1)*(x[9]-(2.0))));
      cinf = 1/(1+exp(-2*(0.4)*(x[9]-(p[0]))));

      taum = cosh ((0.1) * (x[9] - (2.0))/2.0);
      tauc = cosh ((0.4) * (x[9] - (p[0]))/2.0);

      ICa = (4.4) * ninf * (x[9]-(120.0));
      IK = (8.0) * x[10] * (x[9]-(-80.0));
      IKS = (0.15) * x[11] * (x[9]-(-80.0));
      IL = (2.0) * (x[9]-(-60.0));

      aux[9] = (-(ICa + IK + IL + IKS) + (36.0))/(1.2);
      aux[10] = ((4.9) * (minf - x[10])) * taum;
      aux[11] = ((0.005) * (cinf - x[11])) * tauc;
  }
  {
      double ninf, minf, cinf, taum, tauc, ICa, IK, IKS, IL;

      ninf = 1/(1+exp(-2*(0.055)*(x[6]-(-1.2))));
      minf = 1/(1+exp(-2*(0.1)*(x[6]-(2.0))));
      cinf = 1/(1+exp(-2*(0.4)*(x[6]-(p[0]))));

      taum = cosh ((0.1) * (x[6] - (2.0))/2.0);
      tauc = cosh ((0.4) * (x[6] - (p[0]))/2.0);

      ICa = (4.4) * ninf * (x[6]-(120.0));
      IK = (8.0) * x[7] * (x[6]-(-80.0));
      IKS = (0.15) * x[8] * (x[6]-(-80.0));
      IL = (2.0) * (x[6]-(-60.0));

      aux[6] = (-(ICa + IK + IL + IKS) + (36.0))/(1.2);
      aux[7] = ((4.9) * (minf - x[7])) * taum;
      aux[8] = ((0.005) * (cinf - x[8])) * tauc;
  }
  {
      double ninf, minf, cinf, taum, tauc, ICa, IK, IKS, IL;

      ninf = 1/(1+exp(-2*(0.055)*(x[12]-(-1.2))));
      minf = 1/(1+exp(-2*(0.1)*(x[12]-(2.0))));
      cinf = 1/(1+exp(-2*(0.4)*(x[12]-(p[0]))));

      taum = cosh ((0.1) * (x[12] - (2.0))/2.0);
      tauc = cosh ((0.4) * (x[12] - (p[0]))/2.0);

      ICa = (4.4) * ninf * (x[12]-(120.0));
      IK = (8.0) * x[13] * (x[12]-(-80.0));
      IKS = (0.15) * x[14] * (x[12]-(-80.0));
      IL = (2.0) * (x[12]-(-60.0));

      aux[12] = (-(ICa + IK + IL + IKS) + (36.0))/(1.2);
      aux[13] = ((4.9) * (minf - x[13])) * taum;
      aux[14] = ((0.005) * (cinf - x[14])) * tauc;
  }
  aux[0] += -(p[1]) * x[19] * (x[0] - (-70.0)) / (1.2);
  {
     double GV = (0.002) / (1 + exp(-(0.22) * (x[3] - (2.0))));
     aux[19] = (5000.0) * GV * (1-x[19]) - (0.18)*x[19];
  }
  aux[12] += -(p[1]) * x[19] * (x[12] - (-70.0)) / (1.2);
  aux[6] += -(p[1]) * x[19] * (x[6] - (-70.0)) / (1.2);
  aux[6] += -(p[1]) * x[23] * (x[6] - (-70.0)) / (1.2);
  {
     double GV = (0.002) / (1 + exp(-(0.22) * (x[15] - (2.0))));
     aux[23] = (5000.0) * GV * (1-x[23]) - (0.18)*x[23];
  }
  aux[12] += -(p[1]/2.0) * x[23] * (x[12] - (-70.0)) / (1.2);
  aux[3] += -(p[1]/2.0) * x[18] * (x[3] - (-70.0)) / (1.2);
  {
     double GV = (0.002) / (1 + exp(-(0.22) * (x[0] - (2.0))));
     aux[18] = (5000.0) * GV * (1-x[18]) - (0.18)*x[18];
  }
  aux[9] += -(p[1]) * x[18] * (x[9] - (-70.0)) / (1.2);
  aux[0] += -(p[1]) * x[21] * (x[0] - (-70.0)) / (1.2);
  {
     double GV = (0.002) / (1 + exp(-(0.22) * (x[9] - (2.0))));
     aux[21] = (5000.0) * GV * (1-x[21]) - (0.18)*x[21];
  }
  aux[12] += -(p[1]/2.0) * x[21] * (x[12] - (-70.0)) / (1.2);
  aux[15] += -(p[1]) * x[20] * (x[15] - (-70.0)) / (1.2);
  {
     double GV = (0.002) / (1 + exp(-(0.22) * (x[6] - (2.0))));
     aux[20] = (5000.0) * GV * (1-x[20]) - (0.18)*x[20];
  }
  aux[3] += -(p[1]/2.0) * x[20] * (x[3] - (-70.0)) / (1.2);
  aux[3] += -(p[1]) * x[22] * (x[3] - (-70.0)) / (1.2);
  {
     double GV = (0.002) / (1 + exp(-(0.22) * (x[12] - (2.0))));
     aux[22] = (5000.0) * GV * (1-x[22]) - (0.18)*x[22];
  }
  aux[15] += -(p[1]) * x[22] * (x[15] - (-70.0)) / (1.2);
  aux[9] += -(p[1]) * x[22] * (x[9] - (-70.0)) / (1.2);
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
  rop[12] = aux[12];
  rop[13] = aux[13];
  rop[14] = aux[14];
  rop[15] = aux[15];
  rop[16] = aux[16];
  rop[17] = aux[17];
  rop[18] = aux[18];
  rop[19] = aux[19];
  rop[20] = aux[20];
  rop[21] = aux[21];
  rop[22] = aux[22];
  rop[23] = aux[23];
}

void fOLD(double t, double *rop, double *x, double *p, Buffer retard[]) {

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
  aux[0] += -(p[2]) * retard[0](t-p[0]) * (x[0] - (-70.0)) / (1.2);
  aux[3] += -(p[2]) * retard[1](t-p[0]) * (x[3] - (-70.0)) / (1.2);
  aux[6] += -(p[2]) * retard[2](t-p[0]) * (x[6] - (-70.0)) / (1.2);

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

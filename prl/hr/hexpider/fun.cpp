#include <iostream>
#include <cmath>



// SECOND MEMBER FOR RETARDED COPULED NEURONS
void f(double t, double *rop, double *x, double *p) {


    int nvar=24;
    double aux[nvar];


    // p[0] -> b
    // p[1] -> I
    // P[2] -> g_syn
    const double eps = 0.01;
    const double x0 = -1.6;

    // NEURON K
    for(int k=0; k<6; ++k) {
    aux[3*k+0] = x[3*k+1] + x[3*k+0]*x[3*k+0]*(p[0]-x[3*k+0]) - x[3*k+2] + p[1];
    aux[3*k+1] = 1.0 - 5.0*x[3*k+0]*x[3*k+0] - x[3*k+1];
    aux[3*k+2] = eps*(4.0*(x[3*k+0]-x0) - x[3*k+2]);
    }

    // DIFFERENTIAL CONNECTION
    {
      double GV = (0.002) / (1 + exp(-(0.22) * (x[0] - (2.0))));
      aux[18] = (5000.0) * GV * (1-x[18]) - (0.18)*x[18];
    }
    {
      double GV = (0.002) / (1 + exp(-(0.22) * (x[3] - (2.0))));
      aux[19] = (5000.0) * GV * (1-x[19]) - (0.18)*x[19];
    }
    {
      double GV = (0.002) / (1 + exp(-(0.22) * (x[6] - (2.0))));
      aux[20] = (5000.0) * GV * (1-x[20]) - (0.18)*x[20];
    }
    {
      double GV = (0.002) / (1 + exp(-(0.22) * (x[9] - (2.0))));
      aux[21] = (5000.0) * GV * (1-x[21]) - (0.18)*x[21];
    }
    {
      double GV = (0.002) / (1 + exp(-(0.22) * (x[12] - (2.0))));
      aux[22] = (5000.0) * GV * (1-x[22]) - (0.18)*x[22];
    }
    {
      double GV = (0.002) / (1 + exp(-(0.22) * (x[15] - (2.0))));
      aux[23] = (5000.0) * GV * (1-x[23]) - (0.18)*x[23];
    }

    // COUPLING
    // 0 - 3
    // |   |
    // 1 - 4
    // |   |
    // 2 - 5
    
    aux[0] += -(p[2]) * x[19] * (x[0] - (-70.0)) / (1.2);
    aux[0] += -(p[2]) * x[21] * (x[0] - (-70.0)) / (1.2);
    aux[3] += -(p[2]) * x[18] * (x[3] - (-70.0)) / (1.2);
    aux[3] += -(p[2]) * x[20] * (x[3] - (-70.0)) / (1.2);
    aux[3] += -(p[2]) * x[22] * (x[3] - (-70.0)) / (1.2);
    aux[6] += -(p[2]) * x[19] * (x[6] - (-70.0)) / (1.2);
    aux[6] += -(p[2]) * x[23] * (x[6] - (-70.0)) / (1.2);
    aux[9] += -(p[2]) * x[18] * (x[9] - (-70.0)) / (1.2);
    aux[9] += -(p[2]) * x[22] * (x[9] - (-70.0)) / (1.2);
    aux[12] += -(p[2]) * x[19] * (x[12] - (-70.0)) / (1.2);
    aux[12] += -(p[2]) * x[23] * (x[12] - (-70.0)) / (1.2);
    aux[12] += -(p[2]) * x[21] * (x[12] - (-70.0)) / (1.2);
    aux[15] += -(p[2]) * x[20] * (x[15] - (-70.0)) / (1.2);
    aux[15] += -(p[2]) * x[22] * (x[15] - (-70.0)) / (1.2);

    // LOAD FINAL RESULT INTO rop
    for(int k=0; k<nvar; ++k)
      rop[k] = aux[k];
}

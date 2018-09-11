#include <iostream>
#include <cmath>



double sigmoid(double x) {
  return 1.0/(1.0 + exp(80*(x)));
}

// SECOND MEMBER FOR RETARDED COPULED NEURONS
void f(double t, double *rop, double *x, double *p) {

  // p[0] -> vshift
  // p[1] -> I
  // p[2] -> g_syn

  int nvar=18;

  double aux[nvar];

  const double th_syn = -0.03;
  const double E_syn = -0.0625;
  const double E_Na = 0.045;
  const double E_K2 = -0.07;
  const double E_L = -0.046;
  const double g_Na = 160.0;
  const double g_K2 = 30.0;
  const double g_L = 8.0;
  const double C = 0.5;
  const double tau_Na = 0.0405;
  const double tau_K2 = 0.9;

  double I_Na, I_K2, I_L, m_Na;
  double h_Na_inf, m_Na_inf, m_K2_inf;

  for(int i=0; i<6; ++i) {
    h_Na_inf = 1.0/(1.0+exp(500.0*(x[3*i]+0.0325)));
    m_Na_inf = 1.0/(1.0+exp(-150.0*(x[3*i]+0.0305)));
    m_K2_inf = 1.0/(1.0+exp(-83.0*(x[3*i]+0.018+p[0])));

    I_Na = g_Na*m_Na_inf*m_Na_inf*m_Na_inf*x[3*i+11]*(x[3*i]-E_Na);
    I_K2 = g_K2*x[3*i+2]*x[3*i+2]*(x[3*i]-E_K2);
    I_L = g_L*(x[3*i]-E_L);

    aux[3*i] = (-I_Na-I_K2-I_L-p[1])/C;
    aux[3*i+1] = (h_Na_inf-x[3*i+11])/tau_Na;
    aux[3*i+2] = (m_K2_inf-x[3*i+2])/tau_K2;
  }


  // COUPLING


    aux[ 0] += p[2]*(x[ 0]-E_syn) * (sigmoid(x[3]-th_syn) + sigmoid(x[9]-th_syn));
    aux[ 3] += p[2]*(x[ 3]-E_syn) * (0.5*sigmoid(x[0]-th_syn)
                    + 0.5*sigmoid(x[6]-th_syn) + sigmoid(x[12]-th_syn));
    aux[ 6] += p[2]*(x[ 6]-E_syn) * (sigmoid(x[3]-th_syn) + sigmoid(x[15]-th_syn));
    aux[ 9] += p[2]*(x[ 9]-E_syn) * (sigmoid(x[0]-th_syn) + sigmoid(x[12]-th_syn));
    aux[12] += p[2]*(x[12]-E_syn) * (0.5*sigmoid(x[9]-th_syn) + sigmoid(x[3]-th_syn)
                    + 0.5*sigmoid(x[15]-th_syn));
    aux[15] += p[2]*(x[15]-E_syn) * (sigmoid(x[12]-th_syn) + sigmoid(x[6]-th_syn));

  for(int k=0; k<nvar; ++k)
    rop[k] = aux[k];
}

#include <iostream>
#include <cmath>


// SECOND MEMBER FOR RETARDED COPULED NEURONS

void f(double t, double *rop, double *x, double *p) {
  int nvar=3;
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

  // p[0] -> vshift
  // p[1] -> I
  // p[2] -> g_syn

  double I_Na, I_K2, I_L, m_Na;
  double h_Na_inf, m_Na_inf, m_K2_inf;
  // Neurona 0
  // for(int k=0; k<6; ++k) {
    h_Na_inf = 1.0/(1.0+exp(500.0*(x[0]+0.0325)));
    m_Na_inf = 1.0/(1.0+exp(-150.0*(x[0]+0.0305)));
    m_K2_inf = 1.0/(1.0+exp(-83.0*(x[0]+0.018+p[0])));

    I_Na = g_Na*m_Na_inf*m_Na_inf*m_Na_inf*x[1]*(x[0]-E_Na);
    I_K2 = g_K2*x[2]*x[2]*(x[0]-E_K2);
    I_L = g_L*(x[0]-E_L);

    aux[0] = (-I_Na-I_K2-I_L-p[1])/C;
    aux[1] = (h_Na_inf-x[1])/tau_Na;
    aux[2] = (m_K2_inf-x[2])/tau_K2;
  // }


  for(int i=0; i<nvar; ++i) rop[i] = aux[i];

	return;
}

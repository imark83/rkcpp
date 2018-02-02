#include "../variable.hpp"
#include "../rk.hpp"
#include <cstdlib>
#include <cstdio>
#include <cmath>

int nsteps = 0;
int nrejected = 0;

int N = 64;
double yMIN = -0.5;
double yMAX = 0.7;
double YMIN = -0.5;
double YMAX = 0.5;

double H = 1.0/8.0;

double getFLI (const Variable *x) {
  double fli;
  double r, v;
  r = sqrt (x[0][0]*x[0][0] + x[1][0]*x[1][0]);
  v = sqrt (x[2][0]*x[2][0] + x[3][0]*x[3][0]);

  double vv[4];
  for (int i=0; i<4; ++i)
    vv[i] = log ((x[i][1]*x[i][1] + x[i][2]*x[i][2])/r + (x[i][3]*x[i][3] + x[i][4]*x[i][4])/v)/2.0;

  fli = vv[0];
  fli = MAX (fli,vv[1]);
  fli = MAX (fli,vv[2]);
  fli = MAX (fli,vv[3]);

  return fli;
}

int main(int argc, char const *argv[]) {
  std::cout.precision(15);

  int nvar = 4;
  double *fli = (double *) malloc(N*N * sizeof(double));


#pragma omp parallel for
  for(int i=0; i<N; ++i) {
    for(int j=0; j<N; ++j) {
      std::cout << "(" << i << ", " << j << ")" << '\n';
      Variable x[nvar];
      for (int k=0; k<nvar; ++k) x[k] = Variable(4);

    	x[1] = yMIN + ((yMAX-yMIN)*i) / (N-1.0);
    	x[3] = YMIN + ((YMAX-YMIN)*j) / (N-1.0);
    	x[2] = 2.0*H - x[3]*x[3] - x[1]*x[1] + 2.0/3.0*x[1]*x[1]*x[1];

      if (x[2][0] < 0.0) fli[N*j+i] = -5.0;
      else {
        x[2][0] = sqrt (x[2][0]);
      	x[0][1] = x[1][2] = x[2][3] = x[3][4] = 1.0;
        rk(nvar, x, 0, 300, 1, NULL, 1.0e-8);
      	fli[N*j+i] = getFLI(x);
      }
    }
	}


  FILE *fout = fopen("fli.bin", "wb");
  fwrite(fli, sizeof(double), N*N, fout);
  fclose(fout);

  return 0;
}

#include "work.hpp"
#include "rk.hpp"


void work(Task &task) {
  task.result.sn=1000+task.index;
  task.result.period = 1000 - task.vthKS;
  task.result.dutyCycle = 1000 + task.Iext;


  int nvar = 3;
  double x[] = {1.0, 0.5, 0.5};
  double pars[] = {task.vthKS, task.Iext};
  rk(nvar, x, 0, 1000.0, 1000.0, pars, 1.0e-8, 0, task);


}

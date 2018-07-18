#include "work.hpp"
#include "rk.hpp"

#include <iostream>

void work(Task &task, double *x) {
  task.result.sn=0;
  task.result.bn = 0.0;



  int nvar = 27;
  double pars[] = {task.sp, task.km};
  double transient = 2.0 * task.sp;
  double tf = 5.0 * task.sp;
  rk(nvar, x, 0, transient, 100, pars, 1.0e-8, 0, task);
  // rk(nvar, x, 0, 200.0, 200.0, pars, 1.0e-8, 1, task);
  rk(nvar, x, 0, tf, tf, pars, 1.0e-8, 2, task);

}

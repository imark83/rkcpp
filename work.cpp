#include "work.hpp"


void work(Task &task) {
  task.result.sn=1000+task.index;
  task.result.period = 1000 - task.vthKS;
  task.result.dutyCycle = 1000 + task.Iext;
}

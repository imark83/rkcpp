#ifndef __COMMON_HPP__
#define __COMMON_HPP__

#include <stdint.h>

typedef struct {
  uint64_t sn;
  double period;
  double dutyCycle;
  // result_t () : sn(0), period(0.0), dutyCycle(0.0) {}
  // result_t(uint64_t _sn, double _period, double _dutyCycle)
  //       : sn(_sn), period(_period), dutyCycle(_dutyCycle) {}
} result_t;

class Task {
public:
  double vthKS;
  double Iext;
  int index;
  int i;
  int j;
  result_t result;

  Task () {}

  Task(int _index, int _i, int _j, const double *g_vthKS, const double *g_Iext, int M)
      : index(_index), i(_i), j(_j) {
    vthKS = g_vthKS[0];
    if(M>1)
      vthKS += _i*(g_vthKS[1]-g_vthKS[0])/(M-1.0);
    Iext = g_Iext[0];
    if(M>1)
      Iext += _j*(g_Iext[1]-g_Iext[0])/(M-1.0);
  }
};


typedef struct {
  int type;
  int index;
  int chunkSize;
} header_t;


#endif

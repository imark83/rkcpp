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
  double b;
  double I;
  int index;
  int i;
  int j;
  result_t result;

  Task () {}

  Task(int _index, int _i, int _j, const double *g_b, const double *g_I, int M)
      : index(_index), i(_i), j(_j) {
    b = g_b[0];
    if(M>1)
      b += _j*(g_b[1]-g_b[0])/(M-1.0);
    I = g_I[0];
    if(M>1)
      I += _i*(g_I[1]-g_I[0])/(M-1.0);
  }
};


typedef struct {
  int type;
  int index;
  int chunkSize;
} header_t;


#endif

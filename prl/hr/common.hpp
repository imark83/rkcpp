#ifndef __COMMON_HPP__
#define __COMMON_HPP__

#include <stdint.h>

typedef struct {
  uint64_t sn;
  uint64_t bn;
  double orbit[27];
} result_t;

class Task {
public:
  double km;
  double sp;
  int index;
  int i;
  int j;
  result_t result;

  Task () {}

  Task(int _index, int _i, int _j, const double *g_km, const double *g_sp, int M)
      : index(_index), i(_i), j(_j) {
    km = g_km[0];
    if(M>1)
      km += _j*(g_km[1]-g_km[0])/(M-1.0);
    sp = g_sp[0];
    if(M>1)
      sp += _i*(g_sp[1]-g_sp[0])/(M-1.0);
  }
};


typedef struct {
  int type;
  int index;
  int chunkSize;
} header_t;


#endif

#include <iostream>
#include <cstring>
#include <cstdlib>
#include <unistd.h>
#include "buffer.hpp"

Buffer::Buffer()
      : n(0), minPos(0), status(0), defaultV(0.0) {
  T = new double[1];
  V = new double[1];
}

Buffer::Buffer(size_t n, double defaultV)
      : n(n), minPos(0), status(0), defaultV(defaultV) {
  T = new double[n];
  V = new double[n];
}

Buffer::Buffer(const double *T, const double *V, size_t n)
      : n(n), minPos(0), status(n), defaultV(0) {
  this->T = new double[n];
  this->V = new double[n];
  memcpy(this->T, T, n*sizeof(double));
  memcpy(this->V, V, n*sizeof(double));
  minPos = 0;
}

Buffer::Buffer(const Buffer &op)
      : n(op.n), minPos(op.minPos), status(op.status), defaultV(op.defaultV) {
  this->T = new double[n];
  this->V = new double[n];
  memcpy(this->T, op.T, n*sizeof(double));
  memcpy(this->V, op.V, n*sizeof(double));
}


Buffer::~Buffer() {
  delete [] T;
  delete [] V;
}

Buffer & Buffer::operator=(const Buffer &op) {
  this->n = op.n;
  this->minPos = op.minPos;
  this->status = op.status;
  this->defaultV = op.defaultV;
  delete [] this->T;
  delete [] this->V;

  this->T = new double[n];
  this->V = new double[n];
  memcpy(this->T, op.T, n*sizeof(double));
  memcpy(this->V, op.V, n*sizeof(double));
  return (*this);
}

void Buffer::resetTime() {
  double t = T[(minPos-1)%n];
  for(size_t i=0; i<n; ++i) T[i] -= t;
  return;
}



size_t Buffer::getPos(double t) {
  if(t<T[(minPos)] || t>T[(n-1+minPos)%n]) {
    std::cerr << "Value out of Range in Buffer" << std::endl;
    exit(1);
  }
  size_t a=0, b=n-1, c;
  while(b-a > 1) {
    c = (b+a)/2;
    if(T[(c+minPos)%n] < t)
      a = c;
    else
      b = c;
  }
  return a;
}

double Buffer::operator()(double t) {
  if (status < n) return defaultV;
  size_t a = (getPos(t) + minPos)%n;
  size_t b = (a+1)%n;
  // std::cout << a << ", " << b << std::endl;
  // std::cout << V[a] << ", " << V[b] << std::endl;
  double th = (t-T[a])/(T[b]-T[a]);
  return V[a] + th*(V[b]-V[a]);
}
void Buffer::push_back(double t, double v) {
  if (status < n) {
    T[status] = t;
    V[status] = v;
    ++status;
    if(status == n) {
      std::cerr << "buffer filled!" << std::endl;
      // usleep(1000000);
    }
  } else {
    T[minPos] = t;
    V[minPos] = v;
    minPos = (minPos+1)%n;
  }
}

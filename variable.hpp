#ifndef __VARIABLE_HPP__
#define __VARIABLE_HPP__

#include <iostream>
#include <cstdlib>
#include <cmath>

template <class T>
class Variable_ {
public:
  int nder;     // NUMBER OF PARTIAL DERIVATIVES
  T *x;         // VARIABLE AND ITS PARTIAL DERIVATIVES

  // CONSTRUCTORES
  Variable_(int n=0){
    nder = n;
    x = new T[n+1];
    for(int i=0; i<=nder; ++i) x[i] = T(0);
  }
  Variable_(const Variable_ &op){
    nder = op.nder;
    x = new T[nder+1];
    for(int i=0; i<=nder; ++i) x[i] = op.x[i];
  }

  ~Variable_() {
    delete [] x;
  }

  // MEMBER ACCES OPERATORS
  const T &operator[](int i) const {
    if(i>nder){
      std::cerr << "Variable_ index overflow" << std::endl;
      exit(0);
    }
    return x[i];
  }
  T &operator[](int i) {
    if(i>nder){
      std::cerr << "Variable_ index overflow" << std::endl;
      exit(0);
    }
    return x[i];
  }

  // ASSIGNMENT OPERATORS
  Variable_<T> &operator=(const Variable_<T> &op) {
    if (&op == this) return *this;
    if (op.nder != nder){
      nder = op.nder;
      delete [] x;
      x = new T[nder+1];
    }
    for(int i=0; i<nder+1; ++i)
      x[i] = op.x[i];

    return *this;
  }
  Variable_<T> &operator=(T op) {
    x[0] = op;
    for(int i=1; i<nder+1; ++i) x[i] = T(0);

    return *this;
  }


  // NEGATION OPERATORS
  // void neg() {
  //   for(int i=0; i<nder+1; ++i) x[i] = -x[i];
  // }

};


// BASIC OUTPUT
template <class T>
std::ostream &operator<<(std::ostream &output, const Variable_<T> &op) {
  output << op[0];
  if (op.nder) {
    // output << " (" << op[1];
    output << " " << op[1];
    for(int i=2; i<op.nder+1; ++i)
      output << ", " << op[i];
    // output << ")";
  }

  return output;
}



// ARITMETIC OPERATORS

// ADDITION
template <class T>
Variable_<T> operator+(const Variable_<T> &op1, const Variable_<T> &op2) {
  if(op1.nder != op2.nder) {
    std::cerr << "Dimensions missmatch" << std::endl;
    exit(1);
  }
  Variable_<T> rop(op1.nder);
  for(int i=0; i<=rop.nder; ++i)
    rop[i] = op1[i] + op2[i];

  return rop;
}

template <class T>
Variable_<T> operator+(const Variable_<T> &op1, T op2) {
  Variable_<T> rop(op1);
  rop[0] = op1[0] + op2;
  return rop;
}

template <class T>
Variable_<T> operator+(T op2, const Variable_<T> &op1) {
  Variable_<T> rop(op1);
  rop[0] = op1[0] + op2;
  return rop;
}

template <class T>
Variable_<T> operator+(const Variable_<T> &op1, int op2) {
  Variable_<T> rop(op1);
  rop[0] = op1[0] + op2;
  return rop;
}

template <class T>
Variable_<T> operator+(int op2, const Variable_<T> &op1) {
  Variable_<T> rop(op1);
  rop[0] = op1[0] + op2;
  return rop;
}

// SUBTRACTOIN
template <class T>
Variable_<T> operator-(const Variable_<T> &op) {
  Variable_<T> rop(op);
  for(int i=0; i<=op.nder; ++i)
    rop[i] = -rop[i];
  return rop;
}

template <class T>
Variable_<T> operator-(const Variable_<T> &op1, const Variable_<T> &op2) {
  if(op1.nder != op2.nder) {
    std::cerr << "Dimensions missmatch" << std::endl;
    exit(1);
  }
  Variable_<T> rop(op1.nder);
  for(int i=0; i<=rop.nder; ++i)
    rop[i] = op1[i] - op2[i];

  return rop;
}

template <class T>
Variable_<T> operator-(const Variable_<T> &op1, T op2) {
  Variable_<T> rop(op1);
  rop[0] = op1[0] - op2;
  return rop;
}
template <class T>
Variable_<T> operator-(T op2, const Variable_<T> &op1) {
  Variable_<T> rop(op1.nder);
  rop[0] = op2 - op1[0];
  for(int i=1; i<=op1.nder; ++i)
    rop[i] = -op1[i];
  return rop;
}

template <class T>
Variable_<T> operator-(const Variable_<T> &op1, int op2) {
  Variable_<T> rop(op1);
  rop[0] = op1[0] - op2;
  return rop;
}
template <class T>
Variable_<T> operator-(int op2, const Variable_<T> &op1) {
  Variable_<T> rop(op1.nder);
  rop[0] = op2 - op1[0];
  for(int i=1; i<=op1.nder; ++i)
    rop[i] = -op1[i];
  return rop;
}


// PRODUCT
template <class T>
Variable_<T> operator*(const Variable_<T> &op1, const Variable_<T> &op2) {
  if(op1.nder != op2.nder) {
    std::cerr << "Dimensions missmatch" << std::endl;
    exit(1);
  }
  Variable_<T> rop(op1.nder);
  rop[0] = op1[0] * op2[0];
  for (int i=1; i<=rop.nder; ++i)
    rop[i] = op1[0]*op2[i] + op1[i]*op2[0];

  return rop;
}

template <class T>
Variable_<T> operator*(const Variable_<T> &op1, T op2) {
  Variable_<T> rop(op1.nder);
  for(int i=0; i<=rop.nder; ++i)
    rop[i] = op1[i] * op2;
  return rop;
}

template <class T>
Variable_<T> operator*(T op2, const Variable_<T> &op1) {
  Variable_<T> rop(op1.nder);
  for(int i=0; i<=rop.nder; ++i)
    rop[i] = op1[i] * op2;
  return rop;
}

template <class T>
Variable_<T> operator*(const Variable_<T> &op1, int op2) {
  Variable_<T> rop(op1.nder);
  for(int i=0; i<=rop.nder; ++i)
    rop[i] = op1[i] * op2;
  return rop;
}

template <class T>
Variable_<T> operator*(int op2, const Variable_<T> &op1) {
  Variable_<T> rop(op1.nder);
  for(int i=0; i<=rop.nder; ++i)
    rop[i] = op1[i] * op2;
  return rop;
}

template <class T>
Variable_<T> exp(const Variable_<T> &op) {
  Variable_<T> rop(op.nder);
  rop[0] = exp(op[0]);
  for(int i=1; i<=op.nder; ++i)
    rop[i] = exp(op[0]) * op[i];

  return rop;
}

// division
template <class T>
Variable_<T> operator/(const Variable_<T> &op1, const Variable_<T> &op2) {
  if(op1.nder != op2.nder) {
    std::cerr << "Dimensions missmatch" << std::endl;
    exit(1);
  }
  Variable_<T> rop(op1.nder);
  rop[0] = op1[0] / op2[0];
  T aux = op2[0]*op2[0];
  for(int i=1; i<=op1.nder; ++i)
    rop[i] = (op1[i]*op2[0] - op2[i]*op1[0]) / aux;
  return rop;
}

template <class T>
Variable_<T> operator/(T op2, const Variable_<T> &op1) {
  Variable_<T> rop(op1.nder);
  rop[0] = op2/op1[0];
  T aux = op1[0]*op1[0];
  for(int i=1; i<=rop.nder; ++i)
    rop[i] = ((-op2) * op1[i]) / aux;
  return rop;
}

template <class T>
Variable_<T> operator/(const Variable_<T> &op1, T op2) {
  Variable_<T> rop(op1.nder);
  for(int i=0; i<=op1.nder; ++i)
    rop[i] = op1[i]/op2;
  return rop;
}

template <class T>
Variable_<T> operator/(int op2, const Variable_<T> &op1) {
  Variable_<T> rop(op1.nder);
  rop[0] = op2/op1[0];
  T aux = op1[0]*op1[0];
  for(int i=1; i<=rop.nder; ++i)
    rop[i] = ((-op2) * op1[i]) / aux;
  return rop;
}

template <class T>
Variable_<T> operator/(const Variable_<T> &op1, int op2) {
  Variable_<T> rop(op1.nder);
  for(int i=0; i<=op1.nder; ++i)
    rop[i] = op1[i]/op2;
  return rop;
}

template <class T>
Variable_<T> log(const Variable_<T> &op) {
  Variable_<T> rop(op.nder);
  rop[0] = log(op[0]);
  for(int i=1; i<=op.nder; ++i)
    rop[i] = op[i]/op[0];
  return rop;
}

template <class T>
Variable_<T> pow(const Variable_<T> op1, T op2) {
  Variable_<T> rop(op1.nder);
  rop[0] = pow(op1[0],op2);
  for(int i=1; i<=op1.nder; ++i)
    rop[i] = (op2) * pow(op1[0], op2-1.0) * op1[i];

  return rop;
}

template <class T>
Variable_<T> pow(const Variable_<T> op1, int op2) {
  Variable_<T> rop(op1.nder);
  rop[0] = pow(op1[0],op2);
  for(int i=1; i<=op1.nder; ++i)
    rop[i] = (op2) * pow(op1[0], op2-1.0) * op1[i];

  return rop;
}



// FABS (WITH NAN :) )
template <class T>
Variable_<T> fabs(const Variable_<T> &op) {
  Variable_<T> rop(op);
  if(op[0] < 0) for(int i=0; i<=op.nder; ++i)
    rop[i] = -rop[i];
  if(op[0] == 0) for(int i=0; i<=op.nder; ++i)
    rop[i] = T((double) NAN);

  return rop;
}


typedef Variable_<double> Variable;


#endif

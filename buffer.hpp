

class Buffer {
public:
  size_t n;
  double *T;
  double *V;
  size_t minPos;
  size_t status;
  double defaultV;

  Buffer(size_t n, double defaultV);
  Buffer(const double *T, const double *V, size_t n);
  ~Buffer();

  size_t getPos(double t);

  double operator()(double t);
  void push_back(double t, double v);

};

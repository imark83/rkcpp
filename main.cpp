#include <iostream>
#include <fstream>
#include <cmath>
#include <thread>
#include <utility>
#include <deque>
#include <mutex>
#include "rk.hpp"


int nsteps = 0;
int nrejected = 0;

constexpr int M = 10; // Number of points per dimension


using namespace std;


class Task {
public:
    Task(double phi12, double phi13) : phi12{phi12}, phi13{phi13}, result() {}
    Task(Task &&op) {
        phi12 = op.phi12;
        phi13 = op.phi13;
        result = std::move(op.result);
    }

    double phi12;
    double phi13;
    deque<pair<int, double>> result;
};

void worker(double phi12, double phi13, deque<pair<int, double>>& result) {

  int nvar = 24;
  double y[nvar];
  // Mierdas de Buffers

  y[0] = -1.0;
  y[1] = 0.0;
  y[2] = 0.0;
  y[3] = -1.5;
  y[4] = 0.0;
  y[5] = 0.0;
  y[6] = -20.0;
  y[7] = 0.0;
  y[8] = 0.0;
  y[9] = -1.0;
  y[10] = 0.0;
  y[11] = 0.0;
  y[12] = -1.5;
  y[13] = 0.0;
  y[14] = 0.0;
  y[15] = -20.0;
  y[16] = 0.0;
  y[17] = 0.0;

  y[18] = 0.0;
  y[19] = 0.0;
  y[20] = 0.0;
  y[21] = 0.0;
  y[22] = 0.0;
  y[23] = 0.0;


  // START WITH DECOUPLED NETWORK
  double pars[2] = {-26.0, 0.0};
  double poincareThresHold = -30.0;
  rk(nvar, y, 0.0, 4000, 4000,
        pars, 1.0e-8, 0, poincareThresHold);


  double P = rk(nvar, y, 0.0, 1000, 1000,
        pars, 1.0e-8, 1, poincareThresHold);
  std::cerr << "P = " << P << std::endl;

  double x[nvar];     // initial conditions for network
  double z[nvar];     // variables to rock & roll the neuron

  for(int i=0; i<3; ++i) {
    x[9+i] = x[i] = z[i] = y[i];
  }
  x[21] = x[18] = z[18] = y[18];
  rk(nvar, z, 0.0, (P*(1.0-phi12)), (P*(1.0-phi12)),
      pars, 1.0e-8, 0, poincareThresHold);
  for(int i=0; i<3; ++i)
    x[12+i] = x[3+i] = z[i];
  x[22] = x[19] = z[18];


  // GRIND IT... AGAIN!
  for(int i=0; i<3; ++i) {
    z[i] = y[i];
  }
  rk(nvar, z, 0.0, (P*(1.0-phi13)), (P*(1.0-phi13)),
      pars, 1.0e-8, 0, poincareThresHold);

  for(int i=0; i<3; ++i)
    x[15+i] = x[6+i] = z[i];
  x[23] = x[20] = z[18];

  pars[1] = 0.002;
  rk(nvar, x, 0.0, 4000, 4000, pars,
        1.0e-6, 2, poincareThresHold, addressof(result));

  return;
}


int main(int argc, char** argv) {
  vector<Task> tasks;
  tasks.reserve(M*M);

  // Setup tasks
  for(int i = 0; i<M; i++)
    for(int j = 0; j<M; j++)
      tasks.emplace_back(Task(static_cast<double>(i)/M,
              static_cast<double>(j)/M));

  // Initial position for the work
  auto pos = tasks.begin();
  auto end = tasks.end();
  mutex pos_mutex;

  // Setup workers
  vector<thread> workers(thread::hardware_concurrency());
  // vector<thread> workers(1);
  for(auto& th : workers)
    th = thread([&pos, end, &pos_mutex]{
        Task* T;

        while(true) {
          { // Thread safe part
              lock_guard<mutex> lock(pos_mutex);
              if(pos == end) return;
              T = addressof(*pos++);
              //std::cout << std::this_thread::get_id() << " " << T->phi12*M << " " << T->phi13*M << endl;
          }
          // Do the work
          worker(T->phi12, T->phi13, T->result);
        }

      });

  // Wait 'til everything ends
  for(auto& th : workers)
    th.join();

  // OK... print it on screen!
  for(auto &T : tasks) {
    std::cout << T.phi12 << "\t" << T.phi13;
    for(auto& pto : T.result)
      std::cout << "\t" << pto.first << "\t" << pto.second;
    std::cout << std::endl;
  }

  return 0;
}

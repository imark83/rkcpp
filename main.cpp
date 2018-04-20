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

constexpr int M = 30; // Number of points per dimension


using namespace std;


class Task {
public:
    Task(double phi12, double phi13, double retard) : phi12{phi12}, phi13{phi13}, retard{retard}, result() {}
    Task(Task&& op) {
        phi12 = op.phi12;
        phi13 = op.phi13;
        retard = op.retard;
        result = std::move(op.result);
    }

    double phi12;
    double phi13;
    double retard;
    deque<pair<int, double>> result;
};

void worker(double phi12, double phi13, double r, deque<pair<int, double>>& result) {

  int nvar = 12;
  double y[nvar];
  // Mierdas de Buffers
  Buffer retard[] = {Buffer(70000, 0.0), Buffer(70000, 0.0), Buffer(70000, 0.0)};
  Buffer auxBuffer, canonicalBuffer;

  y[0] = -1.0;
  y[1] = 0.0;
  y[2] = 0.0;
  y[3] = -1.5;
  y[4] = 0.0;
  y[5] = 0.0;
  y[6] = -20.0;
  y[7] = 0.0;
  y[8] = 0.0;
  y[9] = 0.0;
  y[10] = 0.0;
  y[11] = 0.0;


  // START WITH DECOUPLED NETWORK
//  double pars[3] = {0, -26.77777777777, 0.0};
  double pars[3] = {10, -26, 0.0};
  double poincareThresHold = -30.0;
  rk(nvar, y, 0.0, 4000, 4000,
        pars, 1.0e-8, 0, poincareThresHold, retard);
  retard[0].resetTime();
  retard[1].resetTime();
  retard[2].resetTime();

  double P = rk(nvar, y, 0.0, 1000, 1000,
        pars, 1.0e-8, 1, poincareThresHold, retard);
  retard[0].resetTime();
  retard[1].resetTime();
  retard[2].resetTime();
//  std::cout << "P = " << P << std::endl;
  canonicalBuffer = retard[0];


  double x[nvar];     // initial conditions for network
  double z[nvar];     // variables to rock & roll the neuron

  for(int i=0; i<3; ++i) {
    x[i] = z[i] = y[i];
  }
  x[9] = z[9] = y[9];
  rk(nvar, z, 0.0, (P*(1.0-phi12)), (P*(1.0-phi12)),
      pars, 1.0e-8, 0, poincareThresHold, retard);
  retard[0].resetTime();
  retard[1].resetTime();
  retard[2].resetTime();
  for(int i=0; i<3; ++i)
    x[3+i] = z[i];
  auxBuffer = retard[0];
  x[10] = z[9];


  // GRIND IT... AGAIN!
  for(int i=0; i<3; ++i) {
    z[i] = y[i];
  }
  retard[0] = canonicalBuffer;
  rk(nvar, z, 0.0, (P*(1.0-phi13)), (P*(1.0-phi13)),
      pars, 1.0e-8, 0, poincareThresHold, retard);
  retard[0].resetTime();
  retard[1].resetTime();
  retard[2].resetTime();

  for(int i=0; i<3; ++i)
    x[6+i] = z[i];
  x[11] = z[9];
  retard[2] = retard[0];
  retard[1] = auxBuffer;
  retard[0] = canonicalBuffer;

  pars[0] = r*P;
  pars[2] = 0.004;
  rk(nvar, x, 0.0, 10000, 10000, pars,
        1.0e-6, 2, poincareThresHold, retard, addressof(result));

  return;
}


int main(int argc, char** argv) {
    vector<Task> tasks;
    tasks.reserve(M*M);

    // Setup tasks
    for(double r = 0.1; r <= 1; r += 0.1)
        for(int i = 0; i<M; i++)
            for(int j = 0; j<M; j++)
                tasks.emplace_back(Task(static_cast<double>(i)/M,
                                        static_cast<double>(j)/M, r));

    // Initial position for the work
    auto pos = tasks.begin();
    auto end = tasks.end();
    mutex pos_mutex;

    // Setup workers
    vector<thread> workers(thread::hardware_concurrency());
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
                worker(T->phi12, T->phi13, T->retard, T->result);
            }

        });

    // Wait 'til everything ends
    for(auto& th : workers)
        th.join();

    // OK... print it on screen!
    for(auto& T : tasks) {
        cout << T.retard << "\t"
             << T.phi12 << "\t"
             << T.phi13;
        for(auto& pto : T.result)
            cout << "\t" << pto.first << "\t" << pto.second;
        cout << endl;
    }



    return 0;
}

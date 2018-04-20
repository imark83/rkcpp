#include <iostream>
#include <fstream>
#include <cmath>
#include <mpi.h>
// #include <utility>
#include <deque>
#include <vector>
#include "rk.hpp"


constexpr int M = 30; // Number of points per dimension


using namespace std;

enum header {TASK_DONE, GIVE_ME, GIVE_YOU, NO_MORE};

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

void work(double phi12, double phi13, deque<pair<int, double>>& result) {

  double r = 0.5;

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


  // SET NEURON (0,0)
  for(int i=0; i<3; ++i) {
    x[i] = z[i] = y[i];
  }
  x[18] = z[18] = y[18];


  // SET NEURON (1,0)
  rk(nvar, z, 0.0, (P*(r)), (P*(r)),
      pars, 1.0e-8, 0, poincareThresHold);
  for(int i=0; i<3; ++i) {
    x[9+i] = z[i];
  }
  z[21] = z[18];



  // SET NEURON (0,1)
  for(int i=0; i<3; ++i){
    z[i] = y[i];
  }
  z[18] = y[18];
  rk(nvar, z, 0.0, (P*(1.0-phi12)), (P*(1.0-phi12)),
      pars, 1.0e-8, 0, poincareThresHold);
  for(int i=0; i<3; ++i)
    x[3+i] = z[i];
  x[19] = z[18];


  // SET NEURON (1,1)
  for(int i=0; i<3; ++i){
    z[i] = y[i];
  }
  z[18] = y[18];
  rk(nvar, z, 0.0, (P*(1.0-phi12+r)), (P*(1.0-phi12+r)),
      pars, 1.0e-8, 0, poincareThresHold);
  for(int i=0; i<3; ++i)
    x[12+i] = z[i];
  x[22] = z[18];



  // SET NEURON (0,2)
  for(int i=0; i<3; ++i) {
    z[i] = y[i];
  }
  z[18] = y[18];
  rk(nvar, z, 0.0, (P*(1.0-phi13)), (P*(1.0-phi13)),
      pars, 1.0e-8, 0, poincareThresHold);
  for(int i=0; i<3; ++i)
    x[6+i] = z[i];
  x[20] = z[18];

  // SET NEURON (1,2)
  for(int i=0; i<3; ++i) {
    z[i] = y[i];
  }
  z[18] = y[18];
  rk(nvar, z, 0.0, (P*(1.0-phi13+r)), (P*(1.0-phi13+r)),
      pars, 1.0e-8, 0, poincareThresHold);
  for(int i=0; i<3; ++i)
    x[15+i] = z[i];
  x[23] = z[18];


  // COUPLE NEURONS AGAIN
  pars[1] = 0.002;
  rk(nvar, x, 0.0, 10000, 10000, pars,
        1.0e-6, 2, poincareThresHold, addressof(result));

  return;
}


int main(int argc, char** argv) {

  // MPI STUFF
  int proc_id, world_size;
  // MPI_Status status;
  int incomingMessage;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  if(world_size < 2) {
    cerr << "need 2 or more slots" << endl;
    exit(1);
  }
  char processor_name[100];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);


  vector<Task> tasks;
  int *buffer = new int[1000];


  if (proc_id == 0) {
    int activeWorkers = world_size-1;
    int pendingTask = taskSize;
    int completedTask = 0;


    tasks.reserve(M*M);
    // Setup tasks
    for(int i = 0; i<M; i++)
      for(int j = 0; j<M; j++)
        tasks.emplace_back(Task(static_cast<double>(i)/M,
              static_cast<double>(j)/M));


    // MAIN COMMUNICATION LOOP
    while(!completedTask || activeWorkers) {
      // QUERY WORKERS ONE BY ONE
      for(int worker=1; worker < world_size; ++worker) {
        MPI_Iprobe(worker, 0, MPI_COMM_WORLD, &incomingMessage,
                &status);
        if(incomingMessage) {
          MPI_Get_count(&status, MPI_INT, &msgLen);
          MPI_Recv(buffer, msgLen, MPI_INT, worker,
                0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        }

      }
    }


  // Initial position for the work
  auto pos = tasks.begin();
  auto end = tasks.end();



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




  MPI_Finalize();

  return 0;
}

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <mpi.h>
#include <utility>
#include <deque>
#include <vector>
#include "rk.hpp"


const int M = 4; // Number of points per dimension


using namespace std;

enum communication_type {TASK_DONE, GIVE_ME, GIVE_YOU, NO_MORE};

class Task {
public:
    Task(double phi21, double phi31) : phi21(phi21), phi31(phi31), result() {}
    Task(const Task &op) :
          phi21(op.phi21), phi31(op.phi31), result(op.result) {}

    double phi21;
    double phi31;
    deque<pair<int, double> > result;
};

typedef struct {
  int type;
  int index;
  double phi21;
  double phi31;
  int ropSize;
} header_t;

void work(double phi21, double phi31, deque<pair<int, double> > &result);

int main(int argc, char** argv) {

  // MPI STUFF
  int proc_id, world_size, bufferLength;
  MPI_Status status;
  int incomingMessage;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  if(world_size < 2) {
    cerr << "need 2 or more slots" << endl;
    MPI_Finalize();
    exit(1);
  }
  char processor_name[100];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);


  vector<Task> tasks;
  char *buffer;



  // ROOT NODE ROUTINES
  if (proc_id == 0) {
    int taskSize = M*M;
    int activeWorkers = world_size-1;
    int pendingTask = taskSize;
    int completedTask = 0;


    tasks.resize(M*M, Task(-1.0, -1.0));
    // Setup tasks
    for(int i = 0; i<M; ++i)
      for(int j = 0; j<M; ++j) {
        Task toEnqueue((1.0*i)/M, (1.0*j)/M);
        tasks[M*i+j] = toEnqueue;
      }


    // MAIN COMMUNICATION LOOP
    while((completedTask != taskSize) || activeWorkers) {
      // cerr << "condition = " << completedTask << endl;
      // QUERY WORKERS ONE BY ONE
      for(int worker=1; worker < world_size; ++worker) {
        MPI_Iprobe(worker, 0, MPI_COMM_WORLD, &incomingMessage,
                &status);
        if(incomingMessage) {
          MPI_Get_count(&status, MPI_BYTE, &bufferLength);
          buffer = new char[bufferLength];
          MPI_Recv(buffer, bufferLength, MPI_CHAR, worker,
                0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          header_t *header = (header_t *) buffer;
          switch (header->type) {
            case TASK_DONE:
              {
              cerr << "proc " << worker << " returns work" << endl;
              ++completedTask;
              vector<pair<int, double> > rop(header->ropSize);
              memcpy(&(rop[0]), buffer+sizeof(header_t),
                  (header->ropSize)*sizeof(pair<int, double>));
              Task toStore(header->phi21, header->phi31);
              for(int i=0; i<rop.size(); ++i)
                toStore.result.push_back(rop[i]);
              tasks[header->index] = toStore;
              delete [] buffer;
              break;
              }
            case GIVE_ME:
              {
              cerr << "proc " << worker << " asks for work" << endl;
              if(pendingTask) {
                int index = taskSize - pendingTask;
                header->type = GIVE_YOU;
                header->index = index;
                header->phi21 = tasks[index].phi21;
                header->phi31 = tasks[index].phi31;
                header->ropSize = 0;
                --pendingTask;
              } else {
                header->type = NO_MORE;
                header->index = -1;
                header->phi21 = -1.0;
                header->phi31 = -1.0;
                header->ropSize = -1;
                --activeWorkers;
              }
              MPI_Send(buffer, sizeof(header_t), MPI_CHAR, worker,
                      0, MPI_COMM_WORLD);
              delete [] buffer;
              break;
              }
          }
        }
      }
    }
  }
  // WORKERS ROUTINES
  else {
    while(1) {
      buffer = new char[sizeof(header_t)];
      header_t *header = (header_t *) buffer;
      header->type = GIVE_ME;
      cerr << "\t\t\tI'm " << processor_name << "-" << proc_id << " and want work" << endl;
      MPI_Send(buffer, sizeof(header_t), MPI_CHAR, 0,
              0, MPI_COMM_WORLD);
      MPI_Recv(buffer, sizeof(header_t), MPI_CHAR, 0,
              0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      if(header->type == NO_MORE) {
        cerr << "\t\t\tI'm " << proc_id << " and recive nothing :( " << buffer[1] << endl;
        cerr << "\t\t\t\t\t\tworker " << processor_name << "-" << proc_id << " finished job" << endl;
        delete [] buffer;
        break;
      }
      cerr << "\t\t\tI'm " << proc_id << " and recive task " << buffer[1] << endl;
      int index = header->index;
      double phi21 = header->phi21;
      double phi31 = header->phi31;
      // DO TASK
      deque<pair<int, double> > result;
      work(phi21, phi31, result);

      vector<pair<int, double> > toSend(result.size());
      for(int i = 0; i < result.size(); ++i)
        toSend[i] = result[i];

      // cerr << "result size = " << result.size() << endl;
      // cerr << "\t\t\t" << proc_id << " gonna send \n\t\t\t";
      // for(int i = 0; i < result.size(); ++i)
      //   cerr << "  (" << toSend[i].first << ", " << toSend[i].second << ")";
      // cerr << endl;

      delete [] buffer;
      buffer = new char[sizeof(header_t)
              + toSend.size()*sizeof(pair<double, int>)];
      header = (header_t *) buffer;

      header->type = TASK_DONE;
      header->index = index;
      header->phi21 = phi21;
      header->phi31 = phi31;
      header->ropSize = toSend.size();
      memcpy(buffer+sizeof(header_t), &(toSend[0]),
              toSend.size()*sizeof(pair<int, double>));
      int bufferLength =
            toSend.size()*sizeof(pair<int, double>) + sizeof(header_t);
      MPI_Send(buffer, bufferLength, MPI_CHAR, 0,
              0, MPI_COMM_WORLD);
      delete [] buffer;
    }
  }


  // ROOT NODE WRITES OUTPUT
  if(proc_id == 0) {
    cerr << "tasks size = " << tasks.size() << endl;
    // OK... print it on screen!
    for(auto &T : tasks) {
      cout << T.phi21 << "\t" << T.phi31;
      for(auto& pto : T.result)
      cout << "\t" << pto.first << "\t" << pto.second;
      cout << endl;
    }
  }

  MPI_Finalize();

  return 0;
}


void work(double phi21, double phi31, deque<pair<int, double> > &result) {

  // result.push_back(make_pair(1, phi21));
  // result.push_back(make_pair(0, phi31));
  // result.push_back(make_pair(2, 2.5));
  //
  // return;

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
  cerr << "P = " << P << endl;

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
  rk(nvar, z, 0.0, (P*(1.0-phi21)), (P*(1.0-phi21)),
      pars, 1.0e-8, 0, poincareThresHold);
  for(int i=0; i<3; ++i)
    x[3+i] = z[i];
  x[19] = z[18];


  // SET NEURON (1,1)
  for(int i=0; i<3; ++i){
    z[i] = y[i];
  }
  z[18] = y[18];
  rk(nvar, z, 0.0, (P*(1.0-phi21+r)), (P*(1.0-phi21+r)),
      pars, 1.0e-8, 0, poincareThresHold);
  for(int i=0; i<3; ++i)
    x[12+i] = z[i];
  x[22] = z[18];



  // SET NEURON (0,2)
  for(int i=0; i<3; ++i) {
    z[i] = y[i];
  }
  z[18] = y[18];
  rk(nvar, z, 0.0, (P*(1.0-phi31)), (P*(1.0-phi31)),
      pars, 1.0e-8, 0, poincareThresHold);
  for(int i=0; i<3; ++i)
    x[6+i] = z[i];
  x[20] = z[18];

  // SET NEURON (1,2)
  for(int i=0; i<3; ++i) {
    z[i] = y[i];
  }
  z[18] = y[18];
  rk(nvar, z, 0.0, (P*(1.0-phi31+r)), (P*(1.0-phi31+r)),
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

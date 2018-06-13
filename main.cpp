#include <iostream>
#include <fstream>
#include <mpi.h>
#include <cstdlib>
#include <vector>


#include "work.hpp"
// #include <cmath>
// #include <cstring>



// #include <deque>
// #include "rk.hpp"



const int M = 2; // Number of points per dimension
const int max_chunkSize = 4; // Max number of allTasks sent to a worker
const double g_vthKS[] = {-28.0, -22};
const double g_Iext[] = {35.0, 36.0};

using namespace std;


enum communication_type {TASK_DONE, GIVE_ME, GIVE_YOU, NO_MORE};




int main(int argc, char** argv) {
  // MPI STUFF
  int proc_id, world_size;
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


  int taskSize = M*M;
  vector<Task> allTasks;
  char *buffer;
  header_t *header;
  Task *task;
  int bufferLen = 0;

  if(proc_id > taskSize) {
    cerr << "\t\t\t NO TASK FOR proc "
          << processor_name << "-" << proc_id << endl;
    MPI_Finalize();
    return 0;
  }

  // ROOT NODE ROUTINES
  if (proc_id == 0) {
    // int taskSize = M*M;
    int activeWorkers = world_size-1;
      if(activeWorkers > taskSize) activeWorkers = taskSize;
    int pendingTask = taskSize;
    int completedTask = 0;


    allTasks.resize(M*M);
    // Setup allTasks
    for(int i = 0; i<M; ++i)
      for(int j = 0; j<M; ++j)
        allTasks[M*i+j] = Task(M*i+j, i, j, g_vthKS, g_Iext, M);


    // MAIN COMMUNICATION LOOP
    while((completedTask != taskSize) || activeWorkers) {
      // cerr << "condition = " << completedTask << endl;
      // QUERY WORKERS ONE BY ONE
      for(int worker=1; worker < world_size; ++worker) {
        MPI_Iprobe(worker, 0, MPI_COMM_WORLD, &incomingMessage,
                &status);
        if(incomingMessage) {
          MPI_Get_count(&status, MPI_BYTE, &bufferLen);
          buffer = new char[bufferLen];
          header = (header_t *) buffer;
          MPI_Recv(buffer, bufferLen, MPI_CHAR, worker,
                0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          switch (header->type) {
            case GIVE_ME: {
              cerr << "proc " << worker << " asks for work" << endl;
              if(pendingTask) {
                cerr << "I have " << pendingTask << " pending tasks" << endl;
                int chunkSize = min(max_chunkSize, pendingTask);
                cerr << "I send you " << chunkSize << " tasks" << endl;

                bufferLen = sizeof(header_t) + chunkSize*sizeof(Task);
                buffer = new char[bufferLen];
                header = (header_t *) buffer;
                task = (Task *) (buffer + sizeof(header_t));

                header->type = GIVE_YOU;
                header->index = taskSize - pendingTask;
                header->chunkSize = chunkSize;
                for(int i=0; i<chunkSize; ++i) {
                  task[i] = allTasks[header->index+i];
                }
                pendingTask -= chunkSize;
              } else {
                header->type = NO_MORE;
                --activeWorkers;
              }

              MPI_Send(buffer, bufferLen, MPI_CHAR, worker,
                    0, MPI_COMM_WORLD);
              delete [] buffer;
              break;
            }
            case TASK_DONE: {
              task = (Task *) (buffer + sizeof(header_t));
              cerr << "proc " << worker << " returns work" << endl;
              completedTask += header->chunkSize;

              for(int i=0; i<header->chunkSize; ++i)
                allTasks[header->index+i].result = task[i].result;

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
      header = (header_t *) buffer;
      header->type = GIVE_ME;
      header->index = -1; header->chunkSize = 0;
      header->chunkSize = 0;
      cerr << "\t\t\tI'm " << processor_name << "-" << proc_id << " and want work" << endl;
      MPI_Send(buffer, sizeof(header_t), MPI_CHAR, 0,
              0, MPI_COMM_WORLD);
      delete [] buffer;

      MPI_Probe(0, 0, MPI_COMM_WORLD, &status);
      MPI_Get_count(&status, MPI_BYTE, &bufferLen);
      buffer = new char[bufferLen];
      header = (header_t *) buffer;
      MPI_Recv(buffer, bufferLen, MPI_CHAR, 0,
              0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      if(header->type == NO_MORE) {
        cerr << "\t\t\tI'm " << processor_name << "-" << proc_id << " and recive nothing :( " << buffer[1] << endl;
        cerr << "\t\t\t\t\t\tworker " << processor_name << "-" << proc_id << " finished job" << endl;
        delete [] buffer;
        break;
      }

      task = (Task *) (buffer + sizeof(header_t));
      cerr << "\t\t\tI'm " << processor_name << "-" << proc_id << " and recive " << header->chunkSize << " indexed at " << header->index << endl;


      // DO TASK
      for(int i=0; i<header->chunkSize; ++i) {
        work(task[i]);
      }

      header->type = TASK_DONE;
      MPI_Send(buffer, bufferLen, MPI_CHAR, 0,
              0, MPI_COMM_WORLD);
      delete [] buffer;
    }
  }


  // ROOT NODE WRITES OUTPUT
  if(proc_id == 0) {
    cerr << "allTasks size = " << allTasks.size() << endl;
    // OK... print it on screen!
    for(int i=0; i<allTasks.size(); ++i) {
      cout << "t[" << i << "] = " << allTasks[i].vthKS << ", " << allTasks[i].Iext << ", " << allTasks[i].result.sn << ", " << allTasks[i].result.period << ", " << allTasks[i].result.dutyCycle << endl;
      // for(auto& pto : T.result)
      // cout << "\t" << pto.first << "\t" << pto.second;
    }
  }

  MPI_Finalize();

  return 0;
}


// void work(double phi12, double phi13, deque<pair<int, double>>& result) {
//
//   double r = 0.5;
//   int nvar = 12;
//   double y[nvar];
//   // Mierdas de Buffers
//   Buffer retard[] = {Buffer(bufferLen, 0.0), Buffer(bufferLen, 0.0), Buffer(bufferLen, 0.0)};
//   Buffer auxBuffer, canonicalBuffer;
//
//   y[0] = -1.0;
//   y[1] = 0.0;
//   y[2] = 0.0;
//   y[3] = -1.5;
//   y[4] = 0.0;
//   y[5] = 0.0;
//   y[6] = -20.0;
//   y[7] = 0.0;
//   y[8] = 0.0;
//   y[9] = 0.0;
//   y[10] = 0.0;
//   y[11] = 0.0;
//
//
//   // START WITH DECOUPLED NETWORK
//   //  double pars[3] = {0, -26.77777777777, 0.0};
//   double pars[3] = {1, vthKS, 0.0};
//   double poincareThresHold = -30.0;
//   rk(nvar, y, 0.0, 4000, 4000,
//         pars, 1.0e-8, 0, poincareThresHold, retard);
//   retard[0].resetTime();
//   retard[1].resetTime();
//   retard[2].resetTime();
//
//   double P = rk(nvar, y, 0.0, 1000, 1000,
//         pars, 1.0e-8, 1, poincareThresHold, retard);
//   retard[0].resetTime();
//   retard[1].resetTime();
//   retard[2].resetTime();
//    std::cerr << "P = " << P << std::endl;
//   canonicalBuffer = retard[0];
//
//
//   double x[nvar];     // initial conditions for network
//   double z[nvar];     // variables to rock & roll the neuron
//
//   for(int i=0; i<12; ++i) {
//     x[i] = z[i] = y[i];
//   }
//   x[9] = z[9] = y[9];
//   rk(nvar, z, 0.0, (P*(1.0-phi12)), (P*(1.0-phi12)),
//       pars, 1.0e-8, 0, poincareThresHold, retard);
//   retard[0].resetTime();
//   retard[1].resetTime();
//   retard[2].resetTime();
//   for(int i=0; i<3; ++i)
//     x[3+i] = z[i];
//   auxBuffer = retard[0];
//   x[10] = z[9];
//
//
//   // GRIND IT... AGAIN!
//   for(int i=0; i<3; ++i) {
//     z[i] = y[i];
//   }
//   retard[0] = canonicalBuffer;
//   rk(nvar, z, 0.0, (P*(1.0-phi13)), (P*(1.0-phi13)),
//       pars, 1.0e-8, 0, poincareThresHold, retard);
//   retard[0].resetTime();
//   retard[1].resetTime();
//   retard[2].resetTime();
//
//   for(int i=0; i<3; ++i)
//     x[6+i] = z[i];
//   x[11] = z[9];
//   retard[2] = retard[0];
//   retard[1] = auxBuffer;
//   retard[0] = canonicalBuffer;
//
//   pars[0] = r*P;
//   pars[2] = 0.004;
//   rk(nvar, x, 0.0, 10000, 10000, pars,
//         1.0e-6, 2, poincareThresHold, retard, addressof(result));
//
//   return;
// }

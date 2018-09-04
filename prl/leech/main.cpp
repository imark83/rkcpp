#include <iostream>
#include <fstream>
#include <mpi.h>
#include <cstdlib>
#include <vector>
#include <unistd.h>


#include "work.hpp"
// #include <cmath>
// #include <cstring>



// #include <deque>
// #include "rk.hpp"



const int M = 64; // Number of points per dimension
const int max_chunkSize = M; // Max number of allTasks sent to a worker
const double g_vshift[] = {-0.03, 0.015};
const double g_I[] = {-0.05, 0.02};

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
        allTasks[M*i+j] = Task(M*i+j, i, j, g_vshift, g_I, M);


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

      double x[] = {-70.0, 0.0, 0.0};

      // DO TASK
      for(int i=0; i<header->chunkSize; ++i) {
        work(task[i], x);
      }

      header->type = TASK_DONE;
      MPI_Send(buffer, bufferLen, MPI_CHAR, 0,
              0, MPI_COMM_WORLD);
      delete [] buffer;
    }
  }


  // ROOT NODE WRITES OUTPUT
  if(proc_id == 0) {
    usleep(100);
    cerr << "allTasks size = " << allTasks.size() << endl;
    // OK... print it on screen!

    int64_t *sn = new int64_t[M*M];
    double *duty = new double[M*M];
    double *period = new double[M*M];


    for(int i=0; i<(int)allTasks.size(); ++i) {
      sn[i] = allTasks[i].result.sn;
      duty[i] = allTasks[i].result.dutyCycle;
      period[i] = allTasks[i].result.period;
      // cout << "t[" << i << "] = " << allTasks[i].vthKS << ", " << allTasks[i].Iext << ", " << allTasks[i].result.sn << ", " << allTasks[i].result.period << ", " << allTasks[i].result.dutyCycle << endl;
      // for(auto& pto : T.result)
      // cout << "\t" << pto.first << "\t" << pto.second;
    }
    std::ofstream fout;

    fout.open("sn.bin", std::ofstream::binary);
    fout.write((const char*) sn, M*M*sizeof(int64_t));
    fout.close();

    fout.open("duty.bin", std::ofstream::binary);
    fout.write((const char*) duty, M*M*sizeof(double));
    fout.close();

    fout.open("period.bin", std::ofstream::binary);
    fout.write((const char*) period, M*M*sizeof(double));
    fout.close();


    delete [] sn;
    delete [] duty;
    delete [] period;
  }

  MPI_Finalize();

  return 0;
}

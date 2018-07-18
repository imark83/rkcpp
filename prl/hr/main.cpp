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



const int M = 8; // Number of points per dimension
const int max_chunkSize = M; // Max number of allTasks sent to a worker
const double g_km[] = {87.0, 88.0};
const double g_sp[] = {1400.0, 1401.0};

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
        allTasks[M*i+j] = Task(M*i+j, i, j, g_km, g_sp, M);


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

      double x[27];
            x[20] = 0.137483;
            x[12] = 0.130489;
            x[21] = 127.498 ;
            x[22] = 125.711 ;
            x[23] = 0.0046091;
            x[6] = 0.991324;
            x[5] = 1.36521e-6;
            x[7] = 3.32817e-7;
            x[9] = 1.43069e-5;
            x[8] = 5.28378e-7;
            x[10] = 0.00865914;
            x[24] = 12.7657;
            x[25] = 13.2176;
            x[4] = 0.597462;
            x[19] = 10.0799;
            x[0] = -87.4094;
            x[2] = 0.991187;
            x[3] = 0.99421;
            x[1] = 0.00103312;
            x[11] = 0.00677689;
            x[13] = 0.0119339;
            x[14] = 0.0664083;
            x[17] = 0.00358545;
            x[18] = 0.995458 ;
            x[15] = 0.00358575 ;
            x[16] = 0.297391 ;
            x[26] = 0.417681;


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

    FILE *fout = fopen("mierda.txt", "w");

    for(int i=0; i<(int)allTasks.size(); ++i) {
      Task &tsk = allTasks[i];
      fprintf(fout, "%e %e %llu %llu",
          tsk.km, tsk.sp, tsk.result.sn, tsk.result.bn);
      for(int j=0; j<27; ++j)
        fprintf(fout, " %.12le", tsk.result.orbit[j]);
      fprintf(fout, "\n");
    }

    fclose(fout);
  }

  MPI_Finalize();

  return 0;
}

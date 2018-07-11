/*
  Preliminary implementation to manage worker ids for OpenMP
  Author: Jun Shirako <shirako@rice.edu>
 */

#ifndef _HABANERO_WORKER_MANAGER
#define _HABANERO_WORKER_MANAGER

#include <iostream>
#include <atomic>
#include <omp.h>

using namespace std;

#define DEFAULT_NUM_WORKERS 24

class WorkerManager {
  int num_workers;

 public:
  WorkerManager() {
  }

  void init(int n) {
    num_workers = n;
  }

  void init() {
    init(DEFAULT_NUM_WORKERS);
  }

  int get_num_workers() {
    return num_workers;
  }

  int get_current_worker_id() {
    int my_id = omp_get_thread_num();
    assert(my_id < num_workers);  // If failed, please specify proper num_workers.
    return my_id;
  }
};

#endif // _HABANERO_WORKER_MANAGER

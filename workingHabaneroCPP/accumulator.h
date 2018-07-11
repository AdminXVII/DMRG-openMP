/*
  Preliminary accumulator implementation specific to array (vector<double>) sum reduction
  Author: Jun Shirako <shirako@rice.edu>

  Reference paper:
    "Finish Accumulators: a Deterministic Reduction Construct for Dynamic Task Parallelism",
    Jun Shirako, Vincent Cave, Jisheng Zhao, Vivek Sarkar.
    The 4th Workshop on Determinism and Correctness in Parallel Programming (WoDet), March 2013.
 */

#ifndef _HABANERO_ACCUMULATOR
#define _HABANERO_ACCUMULATOR

#include <iostream>
#include <vector>
#include <assert.h>
#include <stdlib.h>
#include "TicketLock.h"
#include "WorkerManager.h"

using namespace std;

#define LAZY_POLICY

class accumulator {
  vector<double> result;

#ifdef LAZY_POLICY
  WorkerManager worker_manager;
  vector<double> *local_storage;
#else
  TicketLock ticket_lock;
  vector<double> intermediate_value;
#endif

 public:
  accumulator() {
#ifdef LAZY_POLICY
    worker_manager.init();
    local_storage = new vector<double>[worker_manager.get_num_workers()];
#else
    ticket_lock.init();
#endif
  }

  void put(const vector<double> &val) {
#ifdef LAZY_POLICY
    int i = worker_manager.get_current_worker_id();

    if (local_storage[i].empty()) {
      local_storage[i] = val;

    } else {
      int size1 = local_storage[i].size();
      int size2 = val.size();
      int up = min(size1, size2);
      for (int j = 0; j < up; j++)
	local_storage[i][j] += val[j];
      if (size2 > size1) {
	for (int j = size1; j < size2; j++)
	  local_storage[i].push_back(val[j]);
      }
    }
#else
    ticket_lock.lock();

    int size1 = intermediate_value.size();
    int size2 = val.size();
    int up = min(size1, size2);
    for (int j = 0; j < up; j++)
      intermediate_value[j] += val[j];
    if (size2 > size1) {
      for (int j = size1; j < size2; j++)
	intermediate_value.push_back(val[j]);
    }

    ticket_lock.unlock();
#endif
  }

  vector<double> get() {
    return result;
  }

  void finish() {
#ifdef LAZY_POLICY
    int n = worker_manager.get_num_workers();
    for (int i = 0; i < n; i++) {
      if (local_storage[i].empty())
	continue;

      int size1 = result.size();
      int size2 = local_storage[i].size();
      int up = min(size1, size2);
      for (int j = 0; j < up; j++)
	result[j] += local_storage[i][j];
      if (size2 > size1) {
	for (int j = size1; j < size2; j++)
	  result.push_back(local_storage[i][j]);
      }
      local_storage[i].clear();
    }
#else
    result = intermediate_value;
    intermediate_value.clear();
#endif
  }
};

#endif // _HABANERO_ACCUMULATOR

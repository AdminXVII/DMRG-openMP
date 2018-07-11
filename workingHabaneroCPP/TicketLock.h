/*
  Simplest implementation of busy-wait ticket lock
  Author: Jun Shirako <shirako@rice.edu>

  Reference paper:
    "Algorithms for Scalable Synchronization on Shared-Memory Multiprocessors",
    John M. Mellor-Crummey and Michael L. Scott.
    ACM Transactions on Computer Systems, 9(1):21–65, February 1991.
 */

#ifndef _HABANERO_TICKET_LOCK
#define _HABANERO_TICKET_LOCK

#include <iostream>
#include <atomic>

using namespace std;

class TicketLock {
  std::atomic_uint ticket;
  volatile int current;

 public:
  TicketLock() {
  }

  void init() {
    ticket.store(0);
    current = 0;
  }

  void lock() {
    int my_ticket = ticket.fetch_add(1, memory_order_relaxed);
    while (my_ticket != current);
  }

  void unlock() {
    current = current + 1;
  }
};

#endif // _HABANERO_TICKET_LOCK

#ifndef SUBDOMAIN_ALG_H
#define SUBDOMAIN_ALG_H

#include "GridInfo.hh"

/**
 * This is an abstract base class which defines an interface for any
 * algorithm that produces a grid for a particular node. All
 * processors are intended to run this algorithm. 
 */
class SubdomainAlg
{
private:
protected:
public:
  SubdomainAlg() {}
  virtual ~SubdomainAlg() {}

  /**
   * Subclasses must override this method and implement an algorithm
   * for domain decomposition. Clients should call alloc_grid() on
   * the returned object to perform memory allocation. Memory
   * deallocation occurs automatically when the object falls out of
   * scope due to the destructor.
   *
   * @param rank the rank of the processor we are finding a grid for. 
   * @param size the total number of processors available to us. 
   * @param grid_info an object containing information about the
   * global grid as determined by parsing the input file. 
   *
   * @return a Grid object (class instance) which has its sizes set
   * but which has not yet allocated any memory.
   */
  virtual GridInfo decompose_domain(int rank, int size, 
                                    GridInfo &info) = 0;
};

#endif // SUBDOMAIN_ALG_H

#ifndef SIMPLE_SD_ALG_H
#define SIMPLE_SD_ALG_H

#include "SubdomainAlg.hh"

/**
 *
 */
class SimpleSDAlg : public SubdomainAlg
{
private:
protected:
public:
  SimpleSDAlg();
  virtual ~SimpleSDAlg();

  /**
   * Implements a simple domain decomposition algorithm which divides
   * the domain into an even number (or possibly 3) of blocks, one of
   * which belongs to each processor. 
   *
   * @param rank the rank of the processor we are finding a grid for. 
   * @param size the total number of processors available to us. 
   * @param grid_info an object containing information about the
   * global grid as determined by parsing the input file. 
   *
   * @return a Grid object (class instance) which has its sizes set
   * but which has not yet allocated any memory.
   */
  GridInfo decompose_domain(int rank, int size, const GridInfo &info);  
};

#endif // SIMPLE_SD_ALG_H

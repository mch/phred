#include "BlockResult.hh"

BlockResult::BlockResult()
  : init_(false)
{
  
}

BlockResult::BlockResult(region_t r)
  : region_(r), init_(false)
{
}

BlockResult::~BlockResult()
{

}

void BlockResult::init(const Grid &grid)
{
  MPI_Datatype temp;
  int sizes[3];
  int subsizes[3];
  int starts[3];

  // Setup (convert to local)

  // Create
  MPI_Type_create_subarray(3, &sizes, &subsizes, &starts, 1, &temp);
  MPI_Type_commit(&temp);
  data_.set_datatype(temp);
  data_.set_num(0);

  init_ = true; 
}

Data &BlockResult::get_result(Grid &grid, unsigned int time_step)
{
  if (init_ && result_time(time_step))
  {
    
    data_.set_num(1);
  } 
  else 
  {
    data_.set_num(0);
  }

  return data_;
}

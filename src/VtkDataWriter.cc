#include "VtkDataWriter.hh"

#ifdef USE_VTK

VtkDataWriter::VtkDataWriter(int rank, int size)
  : DataWriter(rank, size)
{
  writer_.SetDataModeToBinary();
}

VtkDataWriter::~VtkDataWriter()
{

}

VtkDataWriter::init()
{

}

VtkDataWriter::deinit()
{

}

VtkDataWriter::add_variable()
{
  
}

unsigned int VtkDataWriter::write_data(unsigned int time_step, 
                                       Data &data, MPI_Datatype t, 
                                       void *ptr, unsigned int len)
{
  
}

#endif 

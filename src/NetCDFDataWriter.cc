#include "NetCDFDataWriter.hh"

NetCDFDataWriter::NetCDFDataWriter(int rank, int size)
  : DataWriter(rank, size), ncid_(0), omode_(NC_WRITE | NC_SHARE),
    fopen_(false)
{
  
}

NetCDFDataWriter::~NetCDFDataWriter()
{
  deinit();
    
}

void NetCDFDataWriter::init()
{
  if (rank_ != 0 || fopen_)
    return;
  
  int status;

  if (filename_.length() > 0) 
  {
    status = nc_open(filename_.c_str(), omode_, &ncid_);
    if (status != NC_NOERR)
    {
      // REMOVE THIS LINE
      cerr << "NetCDFDataWriter error (REMOVE THIS LINE): " 
           << nc_strerror(status) << endl;

      // Maybe we couldn't open it because it DNE, try creating it. 
      status = nc_create(filename_.c_str(), NC_SHARE, &ncid_);
      if (status != NC_NOERR)
        handle_error(status);

    }
  }
  else
  {
    throw exception(); // Must specify a filename. 
  }
}

void NetCDFDataWriter::deinit()
{
  if (rank_ == 0 && fopen_)
  {
    int status = nc_close(ncid_);
    if (status != NC_NOERR)
      handle_error(status);
  }

}

void NetCDFDataWriter::add_variable(Result &result)
{
  int status, dimid;

  if (rank_ != 0)
    return; 

  if (!fopen_)
    throw exception(); // File must be open

  unsigned int num_dims = result.get_num_dims();
  string var_name = result.get_name();

  if (num_dims == 0)
    throw std::exception(); //"Result must have at least one dimension.");

  if (var_name.length() > 0)
  {
    map<string, vector<int>>::iterator iter = dim_ids_.find(var_name);
    if (dim_ids_.end() == iter)
      throw std::exception(); // Duplicates not allowed
  }

  status = nc_redef(ncid_);
  if (status != NC_NOERR)
    handle_error(status);

  unsigned int *lens = result.get_dim_lengths();
  vector<int> dids;

  for (int i = 0; i < num_dims; i++)
  {
    dimid = get_dim(i, lens[i]);
    dids.push_back(dimid);
  }

  dim_ids_[var_name] = dids;

  status = nc_enddef(ncid_);
  if (status != NC_NOERR)
    handle_error(status);
}

void NetCDFDataWriter::handle_data(unsigned int time_step, Data &data)
{

}

int NetCDFDataWriter::get_dim(int i, int size)
{
  string basename, name;
  int status, dimid, idx;
  size_t len;

  switch (i) 
  {
  case 0:
    basename = "x";
    break;
  case 1:
    basename = "y";
    break;
  case 2:
    basename = "z";
    break;
  default:
    basename = "A";
  }

  status = nc_inq_dimid(ncid_, basename.c_str(), &dimid);

  if (status == NC_NOERR) { // Dimension name exists
    status = nc_inq_dimlen(ncid_, dimid, &len);

    if (status != NC_NOERR)
      handle_error(status);

    if (len == size)
      return dimid;
    else {
      idx = 0;
      while (idx < 100) {
        name = basename + idx;
        status = nc_inq_dimid(ncid_, name.c_str(), &dimid);

        if (status == NC_NOERR) {
          status = nc_inq_dimlen(ncid_, dimid, &len);

          if (status == NC_NOERR) 
            handle_error(status);

          if (len == size) 
            return dimid;        
        } else {
          // Create the dimension
          status = nc_def_dim(ncid_, name.c_str(), size, &dimid);
          if (status != NC_NOERR)
            handle_error(status);

          break;
        }
      }
    }
  }
  // Create the dimension
  status = nc_def_dim(ncid_, name.c_str(), size, &dimid);
  if (status != NC_NOERR)
    handle_error(status);

  return dimid;
}

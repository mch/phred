#include "NetCDFDataWriter.hh"
#include <sstream>

using namespace std;

#ifdef USE_NETCDF

NetCDFDataWriter::NetCDFDataWriter(int rank, int size)
  : DataWriter(rank, size), ncid_(0), omode_(NC_WRITE | NC_SHARE),
    fopen_(false)
{}

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

  const vector<int> &dim_lens = result.get_dim_lengths();
  const vector<string> &dim_names = result.get_dim_names();

  string var_name = result.get_name();

  if (dim_lens.size() == 0)
    throw std::exception(); //"Result must have at least one dimension.");

  if (var_name.length() > 0)
  {
    map<string, vector<int> >::iterator iter = dim_ids_.find(var_name);
    if (dim_ids_.end() != iter)
      throw std::exception(); // Duplicates not allowed
  }

  status = nc_redef(ncid_);
  if (status != NC_NOERR)
    handle_error(status);

  // Add variable
  

  // Add dimension id's
  vector<int> dids;
  for (int i = 0; i < dim_lens.size(); i++)
  {
    dimid = get_dim(i, dim_lens[i], dim_names[i]);
    dids.push_back(dimid);
  }

  dim_ids_[var_name] = dids;

  status = nc_enddef(ncid_);
  if (status != NC_NOERR)
    handle_error(status);
}

void *NetCDFDataWriter::write_data(Data &data, MPI_Datatype t, 
                                   void *ptr, unsigned int len)
{
  if (!fopen_)
    throw exception(); // File must be opened and dimensions defined. 

  vector<int> &dids = dim_ids_[data.get_var_name()];
  int status = NC_NOERR;
  int var_id = var_ids_[data.get_var_name()];
  size_t *start, *count;
  
  start = new size_t[dids.size()];
  count = new size_t[dids.size()];

  for (int i = 0; i < dids.size(); i++)
  {
    start[i] = 0;
    count[i] = dids[i]; // This should write the entire region at
                        // once; it's a bit naieve, but it should
                        // work for now...
  }
}

void *NetCDFDataWriter::write_data(int var_id, size_t *start, 
                                   size_t *count, MPI_Datatype t, 
                                   void *ptr, unsigned int len)
{
  int status = NC_NOERR;

  if (t == MPI_CHAR) 
    status = nc_put_vara_schar(ncid_, var_id, start, count, 
                               static_cast<signed char *>(ptr));
  else if (t == MPI_SHORT)
    status = nc_put_vara_short(ncid_, var_id, start, count, 
                               static_cast<short *>(ptr));
  else if (t == MPI_INT)
    status = nc_put_vara_int(ncid_, var_id, start, count, 
                               static_cast<int *>(ptr));
  else if (t == MPI_LONG)
    status = nc_put_vara_long(ncid_, var_id, start, count, 
                               static_cast<long *>(ptr));
  else if (t == MPI_UNSIGNED_CHAR)
    status = nc_put_vara_uchar(ncid_, var_id, start, count, 
                               static_cast<unsigned char *>(ptr));
  else if (t == MPI_UNSIGNED_SHORT) // Not quite right
    status = nc_put_vara_short(ncid_, var_id, start, count, 
                                static_cast<short *>(ptr));
  else if (t == MPI_UNSIGNED) // Not quite right
    status = nc_put_vara_int(ncid_, var_id, start, count, 
                             static_cast<int *>(ptr));
  else if (t == MPI_UNSIGNED_LONG) // Not quite right
    status = nc_put_vara_long(ncid_, var_id, start, count, 
                               static_cast<long *>(ptr));
  else if (t == MPI_FLOAT)
    status = nc_put_vara_float(ncid_, var_id, start, count, 
                               static_cast<float *>(ptr));
  else if (t == MPI_DOUBLE)
    status = nc_put_vara_double(ncid_, var_id, start, count, 
                               static_cast<double *>(ptr));
  else if (t == MPI_LONG_DOUBLE) // Not quite right
    status = nc_put_vara_double(ncid_, var_id, start, count, 
                               static_cast<double *>(ptr));
  //else if (t == MPI_PACKED)
    
  else 
  {
    int num_ints, num_addrs, num_dts, combiner; 

    MPI_Type_get_envelope(t, &num_ints, &num_addrs, &num_dts, 
                          &combiner);

    int *ints;
    MPI_Aint *aints;
    MPI_Datatype *dts;

    ints = new int[num_ints];
    aints = new MPI_Aint[num_addrs];
    dts = new MPI_Datatype[num_dts];

    MPI_Type_get_contents(t, num_ints, num_addrs, num_dts, 
                          ints, aints, dts);

    int i = 0;
    while (i < num_dts && dts[i] > 0)
    {
      ptr = write_data(var_id, start, count, 
                       dts[i], ptr, ints[i]);
      i++;
    }

    delete[] ints;
    delete[] aints;
    delete[] dts;
  }
  
  if (status != NC_NOERR)
    handle_error(status);

  return ptr;
}


int NetCDFDataWriter::get_dim(int i, int size, string basename)
{
  int status, dimid, idx;
  size_t len;
  ostringstream name;

  if (basename.length() == 0)
  {
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
        name.str(basename);
        name <<  idx;
        status = nc_inq_dimid(ncid_, name.str().c_str(), &dimid);

        if (status == NC_NOERR) {
          status = nc_inq_dimlen(ncid_, dimid, &len);

          if (status == NC_NOERR) 
            handle_error(status);

          if (len == size) 
            return dimid;        
        } else {
          // Create the dimension
          status = nc_def_dim(ncid_, name.str().c_str(), size, &dimid);
          if (status != NC_NOERR)
            handle_error(status);

          break;
        }
      }
    }
  } else {
    // Create the dimension
    status = nc_def_dim(ncid_, basename.c_str(), size, &dimid);
    if (status != NC_NOERR)
      handle_error(status);
  }
  return dimid;
}

void NetCDFDataWriter::handle_error(int status)
{
  throw exception();
}

#endif // USE_NETCDF

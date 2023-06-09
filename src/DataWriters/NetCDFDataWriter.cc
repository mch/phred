/* 
   Phred - Phred is a parallel finite difference time domain
   electromagnetics simulator.

   Copyright (C) 2004-2005 Matt Hughes <mhughe@uvic.ca>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
*/

#include "NetCDFDataWriter.hh"

#include <sstream>

#include "../Exceptions.hh"
#include "../Globals.hh"

using namespace std;


#ifdef USE_NETCDF

NetCDFDataWriter::NetCDFDataWriter()
  : ncid_(0), omode_(NC_WRITE | NC_SHARE),
    fopen_(false), clobber_(false)
{}

NetCDFDataWriter::NetCDFDataWriter(const char *filename, 
                                   bool clobber)
  : ncid_(0), omode_(NC_WRITE | NC_SHARE),
    fopen_(false), clobber_(clobber)
{
  set_filename(string(filename));
}


NetCDFDataWriter::~NetCDFDataWriter()
{
  deinit();
}

void NetCDFDataWriter::init()
{
  if (MPI_RANK != rank_ || fopen_)
    return;
  
  int status;

  if (filename_.length() > 0) 
  {
    if (!clobber_)
    {
      status = nc_open(filename_.c_str(), omode_, &ncid_);
      if (status != NC_NOERR)
      {
        // Maybe we couldn't open it because it DNE, try creating it. 
        if (status == 2)
          status = nc_create(filename_.c_str(), NC_SHARE, &ncid_);
        
        if (status != NC_NOERR)
          handle_error(status);
      }
    } 
    else
    {
      status = nc_create(filename_.c_str(), NC_CLOBBER|NC_SHARE, &ncid_);
      
      if (status != NC_NOERR)
        handle_error(status);
    }
    fopen_ = true;
  }
  else
  {
    throw DataWriterException("NetCDFDataWriter: A filename must be "
                              "set before calling init().");
  }
}

void NetCDFDataWriter::deinit()
{
  if (MPI_RANK == rank_ && fopen_)
  {
    int status = nc_close(ncid_);
    if (status != NC_NOERR)
      handle_error(status);

    fopen_ = false;
  }

}

void NetCDFDataWriter::add_variable(Result &result)
{
  int status, dimid, tdim;
  ncdfvar var; 

  if (MPI_RANK != rank_)
    return; 

  if (!fopen_)
    throw DataWriterException("NetCDFDataWriter: Must call init() "
                              "before adding variables.");

  const map<string, Variable *> vars = result.get_variables();
  map<string, Variable *>::const_iterator iter;
  map<string, Variable *>::const_iterator iter_e = vars.end();
  
  for (iter = vars.begin(); iter != iter_e; ++iter)
  {
    Variable *r_var = iter->second;
    const vector<Dimension> &dimensions = r_var->get_dimensions();

    vector<int> dim_lens;
    for(unsigned int idx = 0; idx < dimensions.size(); idx++)
      dim_lens.push_back(dimensions[idx].global_len_);

    var.dim_lens_ = dim_lens;

    var.var_name_ = r_var->get_name();

    if (var.var_name_.size() == 0)
      throw DataWriterException("NetCDFDataWriter: Result variable "
                                "must have a valid "
                                "NetCDF variable name.");
    
    if (var.dim_lens_.size() == 0)
    {
      cerr << "NetCDFDataWriter requires result " << result.get_name()
           << ", variable " << var.var_name_ 
           << " to have at least one dimension.\n";
      throw DataWriterException("NetCDFDataWriter: Variable must have "
                                "at least one dimension.");
    }

    if (var.var_name_.length() > 0)
    {
      map<string, ncdfvar>::iterator iter = vars_.find(var.var_name_);
      if (vars_.end() != iter)
        throw DataWriterException("NetCDFDataWriter: Duplicate variables "
                                  "not allowed");
    }
    
    status = nc_redef(ncid_);
    if (status != NC_NOERR && status != NC_EINDEFINE)
      handle_error(status);
    
    // Add dimension id's
    vector<int> dids;
    for (unsigned int i = 0; i < dimensions.size(); i++)
    {
      string temp = ""; 
      temp = dimensions[i].name_;
      
      dimid = get_dim(i, dimensions[i].global_len_, temp);
      dids.push_back(dimid);
    }
    
    var.time_dim_ = r_var->has_time_dimension();
    if (var.time_dim_)
    {
      status = nc_inq_dimid(ncid_, "t", &tdim);
      
      if (status != NC_NOERR) {
        status = nc_def_dim(ncid_, "t", 
                            NC_UNLIMITED,
                            &tdim);
        
        if (status != NC_NOERR)
          handle_error(status);
      }
      
      dids.insert(dids.begin(), tdim);
      var.dim_lens_.insert(var.dim_lens_.begin(), 0);
    }
    
    var.dim_ids_ = dids;
    
    // Add variable
    int *dimids = new int[dids.size()];
    
    for (unsigned int i = 0; i < dids.size(); i++)
      dimids[i] = dids[i];
    
    // Does the variable already exist? 
    status = nc_inq_varid(ncid_, var.var_name_.c_str(), &var.var_id_);
    if (status == NC_NOERR)
    {
      // Does the variable have the right dimensions? 
      int ndims;
      status == nc_inq_varndims(ncid_, var.var_id_, &ndims);
      if (status != NC_NOERR)
        handle_error(status);
      
      if (static_cast<unsigned int>(ndims) != dids.size())
      {
        cerr << "NetCDFDataWriter found a variable named " 
             << var.var_name_ << " in the file " << filename_
             << ", but it has the wrong number of dimensions. \n";
        throw DataWriterException("NetCDF variable name found, "
                                  "but wrong number of dimensions.");
      }

      status = nc_inq_vardimid(ncid_, var.var_id_, dimids);
      if (status != NC_NOERR)
        handle_error(status);
      
      size_t len = 0;
      unsigned int i = 0;
      
      // Skip the time dimension, because if it exists, it's length
      // won't be zero anyway.
      if (var.time_dim_)
        i++;
      
      for (i = i; i < dids.size(); i++)
      {
        status = nc_inq_dimlen(ncid_, dimids[i], &len);
        
        if (status != NC_NOERR)
          handle_error(status);
        
        if (static_cast<int>(len) != var.dim_lens_[i])
        {
        cerr << "NetCDFDataWriter found a variable named " 
             << var.var_name_ << " in the file " << filename_
             << ", but it has the wrong length for at least one of its "
          "dimensionsions.\n";
          throw DataWriterException("Found a variable with the right "
                                    "name, but one of it's dimensions "
                                    "is the wrong length.");
        }
      }
    } 
    else 
    {
      status = nc_def_var(ncid_, var.var_name_.c_str(), NC_FLOAT, 
                          dids.size(), dimids, &var.var_id_);
    }

    if (status != NC_NOERR)
      handle_error(status);

    delete[] dimids;
    
    status = nc_enddef(ncid_);
    if (status != NC_NOERR)
      handle_error(status);

    // Remember the variable!
    vars_[var.var_name_] = var;
  }

}

unsigned int NetCDFDataWriter::write_data(unsigned int time_step, 
                                          Variable &variable, 
                                          void *ptr, unsigned int len)
{
  if (!fopen_)
    throw DataWriterException("NetCDFDataWriter: File must be opened "
                              "and dimensions defined.");

  const Data &data = variable.get_data();
  string vname = variable.get_name();
  ncdfvar &var = vars_[vname]; // may throw exception

  size_t *start, *count;
  int sz = var.dim_ids_.size();

  start = new size_t[sz];
  count = new size_t[sz];

//   cerr << "NETCDF writer, num dimensions: " << sz << endl;

  for (int i = 0; i < sz; i++)
  {
    start[i] = 0;
    count[i] = var.dim_lens_[i];

//     cerr << "Dim " << i << ", start: " << start[i] << ", count: " 
//          << count[i] << endl;
  }
  
  if (var.time_dim_)
  {
    count[0] = 1;

    map<Variable *, Vardata *>::iterator auxviter 
      = auxvardata_.find(&variable);    
    
    if (auxviter == auxvardata_.end())
      start[0] = 0;
    else
      start[0] = auxviter->second->output_time_;
  }

  unsigned int ret = write_data(var.var_id_, start, count, 
                                variable.get_element_type(), 
                                ptr, len);

  delete[] start;
  delete[] count;

  return ret;
}

unsigned int NetCDFDataWriter::write_data(int var_id, size_t *start, 
                                          size_t *count, MPI_Datatype t, 
                                          void *ptr, unsigned int len)
{
  int status = NC_NOERR;
  unsigned int bytes_written = 0;

  if (t == MPI_CHAR) 
  {
    status = nc_put_vara_schar(ncid_, var_id, start, count, 
                               static_cast<signed char *>(ptr));
    bytes_written = 0;
  }
  else if (t == MPI_BYTE) 
  {
    status = nc_put_vara_uchar(ncid_, var_id, start, count, 
                               static_cast<unsigned char *>(ptr));
    bytes_written = 0;
  }
  else if (t == MPI_SHORT)
  {
    status = nc_put_vara_short(ncid_, var_id, start, count, 
                               static_cast<short *>(ptr));
  }
  else if (t == MPI_INT)
  {
    status = nc_put_vara_int(ncid_, var_id, start, count, 
                               static_cast<int *>(ptr));
  }
  else if (t == MPI_LONG)
  {
    status = nc_put_vara_long(ncid_, var_id, start, count, 
                               static_cast<long *>(ptr));
  }
  else if (t == MPI_UNSIGNED_CHAR)
  {
    status = nc_put_vara_uchar(ncid_, var_id, start, count, 
                               static_cast<unsigned char *>(ptr));
  }
  else if (t == MPI_UNSIGNED_SHORT) // Not quite right
  {
    status = nc_put_vara_short(ncid_, var_id, start, count, 
                                static_cast<short *>(ptr));
  }
  else if (t == MPI_UNSIGNED) // Not quite right
  {
    status = nc_put_vara_int(ncid_, var_id, start, count, 
                             static_cast<int *>(ptr));
  }
  else if (t == MPI_UNSIGNED_LONG) // Not quite right
  {
    status = nc_put_vara_long(ncid_, var_id, start, count, 
                               static_cast<long *>(ptr));
  }
  else if (t == MPI_FLOAT)
  {
    status = nc_put_vara_float(ncid_, var_id, start, count, 
                               static_cast<float *>(ptr));
  }
  else if (t == MPI_DOUBLE)
  {
    status = nc_put_vara_double(ncid_, var_id, start, count, 
                               static_cast<double *>(ptr));
  }
  else if (t == MPI_LONG_DOUBLE) // Not quite right
  {
    status = nc_put_vara_double(ncid_, var_id, start, count, 
                               static_cast<double *>(ptr));
  }
  //else if (t == MPI_PACKED)
    
  else 
  {
//     int num_ints, num_addrs, num_dts, combiner; 

//     MPI_Type_get_envelope(t, &num_ints, &num_addrs, &num_dts, 
//                           &combiner);

//     int *ints;
//     MPI_Aint *aints;
//     MPI_Datatype *dts;

//     ints = new int[num_ints];
//     aints = new MPI_Aint[num_addrs];
//     dts = new MPI_Datatype[num_dts];

//     MPI_Type_get_contents(t, num_ints, num_addrs, num_dts, 
//                           ints, aints, dts);

//     int i = 0;
//     while (i < num_dts && dts[i] > 0)
//     {
//       bytes_written = write_data(var_id, start, count, 
//                                  dts[i], ptr, ints[i]);
//       char *ptr_temp = static_cast<char *>(ptr);
//       ptr_temp += bytes_written;
//       ptr = static_cast<void *>(ptr_temp);
//       i++;
//     }

//     delete[] ints;
//     delete[] aints;
//     delete[] dts;
    throw DataWriterException("NetCDFDataWriter: unknown data type recieved");
  }
  
  if (status != NC_NOERR)
    handle_error(status);

  return bytes_written;
}


int NetCDFDataWriter::get_dim(int i, int size, string basename)
{
  int status, dimid, idx;
  size_t len, sz;

  sz = size; 

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
//     status = nc_inq_dimlen(ncid_, dimid, &len);

//     if (status != NC_NOERR)
//       handle_error(status);

    // Can't allow reuse of dimensions... consider a variable that has
    // 3 dimensions, but which all have the same length. The
    // dimensions must be unique.
 
    //if (len == sz)
    //  return dimid;
    //else 
    {
      idx = 0;
      while (idx < 100) {
        ostringstream name;
        name << basename << idx++;
        status = nc_inq_dimid(ncid_, name.str().c_str(), &dimid);

        if (status == NC_NOERR) {
          status = nc_inq_dimlen(ncid_, dimid, &len);

          if (status != NC_NOERR) 
            handle_error(status);

          //if (len == sz) 
          //  return dimid;        
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
  throw DataWriterException(nc_strerror(status));
}

ostream& NetCDFDataWriter::to_string(ostream &os) const
{
  return os << "NetCDFDataWriter, writing to '"
            << filename_ << "' on rank " << rank_;
}

#endif // USE_NETCDF

#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include <netcdf.h>

using namespace std;

static const unsigned int X_LEN = 2;
static const unsigned int Y_LEN = 8;

void check_status(int status)
{
  if (status != NC_NOERR)
  {
    cout << "A NetCDF error occured: " << status << endl;
    exit(1);
  }
}

// Data varies fastest along x
unsigned int idx(unsigned int x, unsigned int y)
{
  return x + X_LEN * y;
}

int main(int argc, char **argv)
{
  int var_id;
  int ncid;
  int status;
  int xdim, ydim, tdim;
  int dim_ids[3];
  float data[X_LEN * Y_LEN];
  size_t start[] = {0, 0, 0};
  size_t count[] = {1, X_LEN, Y_LEN};

  for (int i = 0; i < X_LEN * Y_LEN; i++)
    data[i] = i * 0.5;

  status = nc_create("test.nc", NC_SHARE, &ncid);
  check_status(status);

  status = nc_def_dim(ncid, "t", NC_UNLIMITED, &tdim);
  check_status(status);

  status = nc_def_dim(ncid, "x", X_LEN, &xdim);
  check_status(status);

  status = nc_def_dim(ncid, "y", Y_LEN, &ydim);
  check_status(status);
  
  dim_ids[0] = tdim;
  dim_ids[1] = xdim;
  dim_ids[2] = ydim;

  status = nc_def_var(ncid, "test", NC_FLOAT, 3, dim_ids, &var_id);
  check_status(status);

  status = nc_enddef(ncid);
  check_status(status);


  for (int k = 0; k < 5; k++)
  {
    for (int i = 0; i < X_LEN * Y_LEN; i++)
      data[i] = i * 0.5 + k;

    status = nc_put_vara_float(ncid, var_id, start, count, data);
    check_status(status);

    start[0]++;
  }



  status = nc_close(ncid);
  check_status(status);

  return 0;
}

#include <iostream>
#include <mpi.h>

using namespace std;

// Calculate the index of a point at coord in a contiguous array
// representing an N dimensional matrix of size size. Data is
// contiguous along the LAST dimension (z in 3 space), as in phred. 
int idx(int N, int *size, int *coord)
{
  int idx = coord[0];

  for (int i = 0; i < N - 1; i++)
    idx = coord[i + 1] + idx * size[i + 1];

  return idx;
}

// Displacement computation recursion helper
int displ_helper(int N, int *dsize, int *lens,
                 int curr_dim, int disp_idx, 
                 int *disps, int *offset,
                 int *negoffset, bool first)
{
  // Base case
  if (curr_dim == N)
  {
    return disp_idx;
  }

  int j = 0;

  if (first) {
    j = 1;

    disp_idx = displ_helper(N, dsize, lens, curr_dim + 1, disp_idx, disps, 
                            offset, negoffset, true);
  }

  for (; j < lens[curr_dim - 1]; j++)
  {
    if (curr_dim == N - 1)
    {
      disps[disp_idx] = disps[disp_idx - *negoffset] + *offset;
      disp_idx++;
    } else {
      disp_idx = displ_helper(N, dsize, lens, curr_dim + 1, disp_idx, disps, 
                              offset, negoffset, false);
    }
  }

  if (first)
  {
    *negoffset = lens[curr_dim - 1];
    *offset *= dsize[curr_dim - 1];
  }

  return disp_idx;
}                          

// Compute displacements for an indexed data type that pulls out a
// chunk of a larger N dimensional matrix. 
void compute_displacements(int N, int *dsize, 
                           int *coord, int *lens,
                           int *num_displacements,
                           int **displacements)
{
  int num_disps = 1;
  int disp_idx = 1;
  int offset = dsize[N - 1];
  int negoffset = 1;

  for (int i = 0; i < N - 1; i++)
  {
    num_disps *= lens[i];
  }

  int *disps = new int [num_disps];

  disps[0] = idx(N, dsize, coord);

  disp_idx = displ_helper(N, dsize, lens, 1, disp_idx, disps, &offset, 
                          &negoffset, true);

  *displacements = disps;
  *num_displacements = num_disps;
}

void idx_test2()
{
  int N = 2;
  int size[] = {5, 6};
  int coord[2];

  for (int i = 0; i < size[0]; i++)
  {
    for (int j = 0; j < size[1]; j++)
    {
      coord[0] = i;
      coord[1] = j;
      cout << "\t" << idx(N, size, coord);
    }
    cout << endl;
  }
  
  // Sub grid extraction...
  coord[0] = 1;
  coord[1] = 1;
  int lens[2] = {2, 2};
  int *displacements;
  int num_disps;

  compute_displacements(N, size, coord, lens, &num_disps, &displacements);

  cout << "Computed displacements, for point (" << coord[0]
       << ", " << coord[1] << "), lens (" << lens[0]
       << ", " << lens[1] << "), got " << num_disps << " back:\n";

  for (int i = 0; i < num_disps; i++)
    cout << "\t" << displacements[i];

  cout << endl << endl;

  delete[] displacements;

  return;
}

void idx_test3()
{
  int N = 3;
  int size[] = {4, 4, 4};
  int coord[3];
  float *matrix = new float[4*4*4];

  for (int i = 0; i < size[2]; i++)
  {
    cout << "K = " << i << endl;
    for (int j = 0; j < size[0]; j++)
    {
      for (int k = 0; k < size[1]; k++)
      {
        coord[2] = i;
        coord[0] = j;
        coord[1] = k;
        
        matrix[idx(N, size, coord)] = idx(N, size, coord);
        cout << "\t" << idx(N, size, coord);
      }
      cout << "\n";
    }
    cout << "--\n";
  }
  
  // Sub grid extraction...
  coord[0] = 1;
  coord[1] = 1;
  coord[2] = 0;
  int lens[3] = {2, 2, 4};
  float *submatrix = new float [2*2*4];
  int *displacements;
  int num_disps;
  MPI_Status status;
  MPI_Datatype submatrix_t;

  compute_displacements(N, size, coord, lens, &num_disps, &displacements);

  cout << "Computed displacements, for point (" << coord[0]
       << ", " << coord[1] << ", " << coord[2] << "), lens (" << lens[0]
       << ", " << lens[1] << ", " << lens[2] << ") got " 
       << num_disps << " back:\n";

  int *lengths = new int [num_disps];
  for (int i = 0; i < num_disps; i++)
  {
    lengths[i] = 4;
    cout << "\t" << displacements[i] << ", " << lengths[i];
  }

  cout << endl << endl;

  MPI_Type_indexed(num_disps, lengths, displacements, MPI_FLOAT, &submatrix_t);
  MPI_Type_commit(&submatrix_t);
  
  MPI_Sendrecv(matrix, 1, submatrix_t, 0, 1, 
               submatrix, 4*2*2, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &status);

  cout << "Submatrix: \n";
  for (int i = 0; i < 4; i++)
  {
    cout << "K = " << i << endl;
    for (int j = 0; j < 2; j++)
    {
      for (int k = 0; k < 2; k++)
      {
        coord[0] = j;
        coord[1] = k;
        coord[2] = i;

        cout << "\t" << idx(N, lens, coord) << "=" 
             << submatrix[idx(N, lens, coord)];
      }
      cout << "\n";
    }
    cout << "--\n";
  }
  
  cout << endl << endl;

  MPI_Type_free(&submatrix_t); 

  delete[] displacements;
  delete[] submatrix;
  delete[] matrix;
  return;
}

int main (int argc, char **argv)
{
  int rank, size;
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
    
  float *matrix = new float[4*4];
  float *submatrix = new float[2*2];

  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      matrix[j + i*4] = (i+1) * (j+2) * 0.5;

  // Use an indexed type to extract the diagonal
  float *diagonal = new float[4];
  int lengths[4] = {1, 1, 1, 1};
  int disps[4] = {0, 5, 10, 15};
  
  MPI_Datatype diag;
  MPI_Type_indexed(4, lengths, disps, MPI_FLOAT, &diag);
  MPI_Type_commit(&diag);
  
  MPI_Sendrecv(matrix, 1, diag, 0, 1, 
               diagonal, 4, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &status);

  cout << "Initial matrix: \n";
  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
      cout << "\t" << matrix[j + i*4];
    cout << endl;
  }

  cout << "\n\ndiagonal: \n";
  for (int i = 0; i < 4; i++)
    cout << "\t" << diagonal[i];
  cout << endl;

  MPI_Datatype submatrix_t;
  lengths[0] = 2;
  lengths[1] = 2;

  int *displacements;
  int num_disps;
  int lens[] = {2, 2};
  int dsize[] = {4, 4};
  int coord[] = {0, 2};

  compute_displacements(2, dsize, coord, lens, &num_disps, &displacements);

  cout << "Computed displacements, for point (" << coord[0]
       << ", " << coord[1] << "), lens (" << lens[0]
       << ", " << lens[1] << "), got " << num_disps << " back:\n";

  for (int i = 0; i < num_disps; i++)
    cout << "\t" << displacements[i];

  cout << endl << endl;

  MPI_Type_indexed(2, lens, displacements, MPI_FLOAT, &submatrix_t);
  MPI_Type_commit(&submatrix_t);
  
  MPI_Sendrecv(matrix, 1, submatrix_t, 0, 1, 
               submatrix, 4, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &status);

  cout << "\nsubmatrix: \n";
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
      cout << "\t" << submatrix[j + i*2];
    cout << endl;
  }
  cout << endl;
  cout << endl;
  
  MPI_Type_free(&submatrix_t); 
  MPI_Type_free(&diag);
  delete[] diagonal;
  delete[] matrix;
  delete[] submatrix;

  cerr << "idx_test2:" << endl;
  idx_test2();

  cerr << "idx_test3:" << endl;
  idx_test3();

  MPI_Finalize();

  return 0;
}

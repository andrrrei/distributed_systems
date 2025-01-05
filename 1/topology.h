#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include <mpi.h>

#define N_ROW 8
#define N_COL 8


MPI_Comm create_2d_topology(int *rank);

void get_2d_coords(MPI_Comm comm, int rank, int coords[2]);

#endif

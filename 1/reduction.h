#ifndef REDUCTION_H
#define REDUCTION_H

#include <mpi.h>

void reduce_2d_rc(MPI_Comm comm, int coords[2], int *data, int *res_rank);

void broadcast_2d_rc(MPI_Comm comm, int coords[2], int *data, int *res_rank);

#endif

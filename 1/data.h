#ifndef DATA_H
#define DATA_H

#include <mpi.h>

int generate_value(int rank);

void gather_and_print_generated_data(MPI_Comm comm, int rank, int coords[2], int data);

void gather_and_print_final_data(MPI_Comm comm,
                                 int rank,
                                 int res_rank,
                                 int res_coords[2],
                                 int best_data);

#endif

#include "topology.h"
#include <stdio.h>

/*
Инициализирует MPI, создаёт двумерный коммуникатор для размерности N_ROW x N_COL.
*/
MPI_Comm create_2d_topology(int *rank) {
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, rank);

    int dims[2] = {N_ROW, N_COL};
    int periods[2] = {0, 0}; // Без закольцовки
    int reorder = 0;         // Не менять ранги автоматически

    MPI_Comm comm_2d;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &comm_2d);

    return comm_2d;
}

/*
Возвращает двумерные координаты процесса в созданной топологии
*/
void get_2d_coords(MPI_Comm comm, int rank, int coords[2]) {
    MPI_Cart_coords(comm, rank, 2, coords);
}

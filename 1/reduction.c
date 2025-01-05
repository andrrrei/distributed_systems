#include "reduction.h"
#include "topology.h"
#include <mpi.h>
#include <stdio.h>

/*
Возвращает rank процесса с координатами (row, col).
*/
static int get_neighbour_rank(MPI_Comm comm, int row, int col) {
    if (row < 0 || row >= N_ROW || col < 0 || col >= N_COL) {
        return -1;
    }
    int neighbour_rank;
    int coords[2] = { row, col };
    MPI_Cart_rank(comm, coords, &neighbour_rank);
    return neighbour_rank;
}

/**
Редукция по строкам, потом по столбцам.
- Сначала каждый лидер строки (col=0) собирает max из этой строки;
- Затем процесс (0,0) собирает max по всем лидерам строк.
*/
void reduce_2d_rc(MPI_Comm comm, int coords[2],
                  int *data, int *res_rank) 
{
    int row = coords[0];
    int col = coords[1];

    // ==============================
    // 1) Сбор max в лидера строки (col=0)
    // ==============================
    if (col == 0) {
        for (int c = 1; c < N_COL; c++) {
            int source_rank = get_neighbour_rank(comm, row, c);
            int msg[2];
            MPI_Recv(msg, 2, MPI_INT, source_rank, 0, comm, MPI_STATUS_IGNORE);
            int recv_res_rank = msg[0];
            int recv_data      = msg[1];

            // Логика "maxloc"
            if (recv_data > *data) {
                *data = recv_data;
                *res_rank = recv_res_rank;
            } else if (recv_data == *data && recv_res_rank < *res_rank) {
                *res_rank = recv_res_rank;
            }
        }
    } else {
        int leader_rank = get_neighbour_rank(comm, row, 0);
        int msg[2] = { *res_rank, *data };
        MPI_Send(msg, 2, MPI_INT, leader_rank, 0, comm);
    }

    // ==============================
    // 2) Сбор max у (0,0) среди лидеров строк
    // ==============================
    if (col == 0) {
        if (row == 0) {
            for (int r = 1; r < N_ROW; r++) {
                int source_rank = get_neighbour_rank(comm, r, 0);
                int msg[2];
                MPI_Recv(msg, 2, MPI_INT, source_rank, 0, comm, MPI_STATUS_IGNORE);
                int recv_res_rank = msg[0];
                int recv_data      = msg[1];

                // Логика "maxloc"
                if (recv_data > *data) {
                    *data = recv_data;
                    *res_rank = recv_res_rank;
                } else if (recv_data == *data && recv_res_rank < *res_rank) {
                    *res_rank = recv_res_rank;
                }
            }
        } else {
            int global_leader = get_neighbour_rank(comm, 0, 0);
            int msg[2] = { *res_rank, *data };
            MPI_Send(msg, 2, MPI_INT, global_leader, 0, comm);
        }
    }
}

/*
Broadcast по столбцам, потом по строкам
*/
void broadcast_2d_rc(MPI_Comm comm, int coords[2],
                     int *data, int *res_rank)
{
    int row = coords[0];
    int col = coords[1];

    // ==============================
    // 1) Рассылка лидерам строк по столбцу (col=0).
    // ==============================
    if (col == 0) {
        if (row == 0) {
            for (int r = 1; r < N_ROW; r++) {
                int target = get_neighbour_rank(comm, r, 0);
                int msg[2] = { *res_rank, *data };
                MPI_Send(msg, 2, MPI_INT, target, 0, comm);
            }
        } else {
            int global_leader = get_neighbour_rank(comm, 0, 0);
            int msg[2];
            MPI_Recv(msg, 2, MPI_INT, global_leader, 0, comm, MPI_STATUS_IGNORE);
            *res_rank = msg[0];
            *data      = msg[1];
        }
    }

    // ==============================
    // 2) Рассылка по строке
    // ==============================
    if (col == 0) {
        for (int c = 1; c < N_COL; c++) {
            int target = get_neighbour_rank(comm, row, c);
            int msg[2] = { *res_rank, *data };
            MPI_Send(msg, 2, MPI_INT, target, 0, comm);
        }
    } else {
        int leader = get_neighbour_rank(comm, row, 0);
        int msg[2];
        MPI_Recv(msg, 2, MPI_INT, leader, 0, comm, MPI_STATUS_IGNORE);
        *res_rank = msg[0];
        *data      = msg[1];
    }
}

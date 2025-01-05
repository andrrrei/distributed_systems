#include "data.h"
#include "topology.h"
#include <stdlib.h>
#include <stdio.h>


int generate_value(int rank) {
    srand(rank);
    return rand() % 1000000; // ограничим диапазон
}

/*
Сбор и печать исходных данных
*/
void gather_and_print_generated_data(MPI_Comm comm, int rank, int coords[2], int data) {
    int size;
    MPI_Comm_size(comm, &size);

    // Собираем 4 числа от каждого процесса: rank, row, col, data
    int localArr[4] = { rank, coords[0], coords[1], data };

    // Массив для приёма — только в rank=0
    int *globalArr = NULL;
    if (rank == 0) {
        globalArr = (int*)malloc(4 * size * sizeof(int));
    }

    // Сбор данных
    MPI_Gather(localArr, 4, MPI_INT,
               globalArr, 4, MPI_INT,
               0, comm);

    // Печатаем таблицу
    if (rank == 0) {
        printf("GENERATED DATA (all processes):\n");
        printf("┌─────┬─────────┬─────────┐\n");
        printf("│ PID │ coords  │ data    │\n");
        printf("├─────┼─────────┼─────────┤\n");

        for (int i = 0; i < size; i++) {
            int offset = i * 4;
            int pRank = globalArr[offset + 0];
            int r     = globalArr[offset + 1];
            int c     = globalArr[offset + 2];
            int d     = globalArr[offset + 3];

            printf("│ %3d │ (%2d,%2d) │ %7d │\n", pRank, r, c, d);
        }

        printf("└─────┴─────────┴─────────┘\n\n");
        free(globalArr);
    }
}

/*
Сбор и печать итоговых данных
*/
void gather_and_print_final_data(MPI_Comm comm,
                                 int rank,
                                 int res_rank,
                                 int res_coords[2],
                                 int best_data) {
    int size;
    MPI_Comm_size(comm, &size);

    int localArr[5] = { rank,
                        res_rank,
                        res_coords[0],
                        res_coords[1],
                        best_data };

    // Массив приёма в rank=0
    int *globalArr = NULL;
    if (rank == 0) {
        globalArr = (int*)malloc(5 * size * sizeof(int));
    }

    MPI_Gather(localArr, 5, MPI_INT,
               globalArr, 5, MPI_INT,
               0, comm);

    if (rank == 0) {
        printf("FINAL RESULTS (each process sees the same global max):\n");
        printf("┌─────┬────────────┬──────────────┬───────────┐\n");
        printf("│ PID │ res_rank   │ res_coords   │ rea_data  │\n");
        printf("├─────┼────────────┼──────────────┼───────────┤\n");

        for (int i = 0; i < size; i++) {
            int offset = i * 5;
            int pRank  = globalArr[offset + 0];
            int bRank  = globalArr[offset + 1];
            int brRow  = globalArr[offset + 2];
            int brCol  = globalArr[offset + 3];
            int bData  = globalArr[offset + 4];
            printf("│ %3d │ %9d  │ (%2d,%2d)      │ %9d │\n",
                   pRank, bRank, brRow, brCol, bData);
        }

        printf("└─────┴────────────┴──────────────┴───────────┘\n\n");
        free(globalArr);
    }
}

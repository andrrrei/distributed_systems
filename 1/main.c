#include <stdio.h>
#include <mpi.h>

#include "topology.h"
#include "data.h"
#include "reduction.h"

int main(int argc, char *argv[]) {
    int rank;
    // Создаём 2D-топологию и инициализируем MPI
    MPI_Comm comm_2d = create_2d_topology(&rank);

    // Узнаем координаты (row, col)
    int coords[2];
    get_2d_coords(comm_2d, rank, coords);

    // Генерируем для каждого процесса случайное число
    int data = generate_value(rank);

    int res_rank = rank;

    gather_and_print_generated_data(comm_2d, rank, coords, data);

    // 1) Найдём глобальный максимум и ранг его источника
    reduce_2d_rc(comm_2d, coords, &data, &res_rank);
    MPI_Barrier(comm_2d);

    // 2) Рассылаем итоговые данные всем
    broadcast_2d_rc(comm_2d, coords, &data, &res_rank);
    MPI_Barrier(comm_2d);

    // Узнаем координаты процесса, чей res_rank
    int res_coords[2];
    get_2d_coords(comm_2d, res_rank, res_coords);

    // Печатаем результат
    gather_and_print_final_data(comm_2d, rank, res_rank, res_coords, data);
    
    MPI_Finalize();
    return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include <mpi-ext.h> 

#include <signal.h>

#define Max(a, b) ((a) > (b) ? (a) : (b))

// ---------------- ПАРАМЕТРЫ ЗАДАЧИ ----------------
#define N (4096 + 2)
double maxeps = 0.1e-7;
int itmax = 100;

// ------------------ ГЛОБАЛЬНЫЕ ПЕРЕМЕННЫЕ ------------------
double A[N][N], B[N][N];
double eps;
int rank_world, size_world;
MPI_Comm main_comm;
MPI_Errhandler errh;

int iteration_to_crash = -1; 
int rank_to_crash = -1;
int it = 1;

static int has_died = 0; // Флаг, что у нас уже был сбой

// ------------------ ПРОТОТИПЫ ФУНКЦИЙ ------------------
void init(int load_from_checkpoint, int *iteration_start);
void relax();
void verify();
void save_checkpoint(int iteration);
int  load_checkpoint(int *iteration);
void error_handler(MPI_Comm *comm, int *error_code, ...);

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_world);
    MPI_Comm_size(MPI_COMM_WORLD, &size_world);

    // Создаём обработчик ошибок
    MPI_Comm_create_errhandler(error_handler, &errh);

    main_comm = MPI_COMM_WORLD;
    // Привязываем обработчик к главному коммуникатору
    MPI_Comm_set_errhandler(main_comm, errh);

    if (argc == 3) {
        iteration_to_crash = atoi(argv[1]);
        rank_to_crash      = atoi(argv[2]);
    }

    int iteration_start = 1;
    // Пытаемся загрузить чекпоинт (если есть)
    init(/*load_from_checkpoint=*/1, &iteration_start);

    // ------------------ Основной цикл итераций ------------------
    while (it <= itmax) {
        eps = 0.0;

        // Имитация сбоя
        if (rank_world == rank_to_crash && it == iteration_to_crash && !has_died) {
            raise(SIGKILL);
        }

        relax();

        if (rank_world == 0) {
            printf("Iteration = %4d\n", it);
        }

        if ((it % 10) == 0) {
            save_checkpoint(it);
        }

        MPI_Barrier(main_comm);

        it++;
    }


    verify();

    MPI_Finalize();
    return 0;
}

// ========================================================================
//                 ФУНКЦИЯ ИНИЦИАЛИЗАЦИИ И ЗАГРУЗКИ
// ========================================================================
void init(int load_from_checkpoint, int *iteration_start)
{
    if (load_from_checkpoint) {
        int loaded = load_checkpoint(iteration_start);
        if (loaded) {
            if (rank_world == 0) {
                printf("Checkpoint loaded. Continue from iteration %d\n", *iteration_start);
            }
            return;
        }
    }
    if (rank_world == 0) {
        printf("No valid checkpoint found. Initialize from scratch.\n");
    }
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (i==0 || i==N-1 || j==0 || j==N-1)
                A[i][j] = 0.0;
            else
                A[i][j] = (1.0 + i + j);
        }
    }
    *iteration_start = 1;
}

// ========================================================================
//           ФУНКЦИЯ ОДНОГО ШАГА РЕЛАКСАЦИИ
// ========================================================================
void relax()
{
    int rank_size;
    MPI_Comm_rank(main_comm, &rank_world);
    MPI_Comm_size(main_comm, &rank_size);

    int start_i =  rank_world    * (N / rank_size);
    int end_i   = (rank_world+1) * (N / rank_size);

    if (rank_world == rank_size - 1) {
        end_i = N;
    }

    // 1) Вычисляем новую матрицу B
    for (int i = start_i + 1; i < end_i - 1; i++)
    {
        for (int j = 1; j < N - 1; j++)
        {
            B[i][j] = 0.25 * ( A[i-1][j] + A[i+1][j]
                             + A[i][j-1] + A[i][j+1] );
        }
    }

    // 2) Подсчёт локального eps и обновление A
    double local_eps = 0.0;
    for (int i = start_i + 1; i < end_i - 1; i++)
    {
        for (int j = 1; j < N - 1; j++)
        {
            double diff = fabs(A[i][j] - B[i][j]);
            A[i][j] = B[i][j];
            local_eps = Max(local_eps, diff);
        }
    }

    // 3) Собираем глобальный eps
    MPI_Allreduce(&local_eps, &eps, 1, MPI_DOUBLE, MPI_MAX, main_comm);

    // 4) Обмен граничными строками
    if (rank_world < rank_size - 1) {
        MPI_Send(&A[end_i-2][0], N, MPI_DOUBLE, rank_world+1, 0, main_comm);
        MPI_Recv(&A[end_i-1][0], N, MPI_DOUBLE, rank_world+1, 0, main_comm, MPI_STATUS_IGNORE);
    }
    if (rank_world > 0) {
        MPI_Recv(&A[start_i][0], N, MPI_DOUBLE, rank_world-1, 0, main_comm, MPI_STATUS_IGNORE);
        MPI_Send(&A[start_i+1][0], N, MPI_DOUBLE, rank_world-1, 0, main_comm);
    }
}

// ========================================================================
//          ФУНКЦИЯ СОХРАНЕНИЯ СОСТОЯНИЯ
// ========================================================================
void save_checkpoint(int iteration)
{
    char fname[128] = "checkpoint_iter_last.dat";

    MPI_File fh;
    int rc = MPI_File_open(main_comm, fname,
                           MPI_MODE_WRONLY | MPI_MODE_CREATE,
                           MPI_INFO_NULL, &fh);
    if (rc != MPI_SUCCESS) {
        if (rank_world == 0) {
            fprintf(stderr, "Error opening checkpoint file: %s\n", fname);
        }
        return;
    }

    MPI_Offset offset = 0;
    if (rank_world == 0) {
        MPI_File_write_at(fh, offset, &iteration, 1, MPI_INT, MPI_STATUS_IGNORE);
    }
    offset += sizeof(int); 

    int rank_size;
    MPI_Comm_size(main_comm, &rank_size);

    int start_i =  rank_world    * (N / rank_size);
    int end_i   = (rank_world+1) * (N / rank_size);
    if (rank_world == rank_size - 1) {
        end_i = N;
    }
    int local_height = end_i - start_i;

    // Вычисляем смещение (сколько элементов пропустить)
    MPI_Offset block_offset_elems = 0;
    for (int r = 0; r < rank_world; r++) {
        int r_start =  r    * (N / rank_size);
        int r_end   = (r+1) * (N / rank_size);
        if (r == rank_size - 1) r_end = N;
        int r_height = r_end - r_start;
        block_offset_elems += (r_height * N);
    }
    MPI_Offset file_offset = offset + block_offset_elems * sizeof(double);

    MPI_File_write_at(fh, file_offset,
                      &A[start_i][0],
                      local_height * N,
                      MPI_DOUBLE,
                      MPI_STATUS_IGNORE);

    MPI_File_close(&fh);

    if (rank_world == 0) {
        printf("Checkpoint saved at iteration %d\n", iteration);
    }
}

// ========================================================================
//          ФУНКЦИЯ ЗАГРУЗКИ СОСТОЯНИЯ
// ========================================================================
int load_checkpoint(int *iteration)
{
    char fname[128] = "checkpoint_iter_last.dat";

    MPI_File fh;
    if (MPI_File_open(main_comm, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh)
        != MPI_SUCCESS) {
        return 0;
    }

    MPI_Offset offset = 0;
    int tmp_it = 0;
    if (rank_world == 0) {
        MPI_File_read_at(fh, offset, &tmp_it, 1, MPI_INT, MPI_STATUS_IGNORE);
    }
    // Рассылаем значение остальным
    MPI_Bcast(&tmp_it, 1, MPI_INT, 0, main_comm);
    offset += sizeof(int);

    int rank_size;
    MPI_Comm_size(main_comm, &rank_size);

    int start_i =  rank_world    * (N / rank_size);
    int end_i   = (rank_world+1) * (N / rank_size);
    if (rank_world == rank_size - 1) {
        end_i = N;
    }
    int local_height = end_i - start_i;

    MPI_Offset block_offset_elems = 0;
    for (int r = 0; r < rank_world; r++) {
        int r_start =  r    * (N / rank_size);
        int r_end   = (r+1) * (N / rank_size);
        if (r == rank_size - 1) r_end = N;
        int r_height = r_end - r_start;
        block_offset_elems += (r_height * N);
    }

    MPI_Offset file_offset = offset + block_offset_elems * sizeof(double);

    MPI_File_read_at(fh, file_offset,
                     &A[start_i][0],
                     local_height * N,
                     MPI_DOUBLE,
                     MPI_STATUS_IGNORE);

    MPI_File_close(&fh);

    // Возвращаем номер итерации
    *iteration = tmp_it + 1;

    return 1;
}

// ========================================================================
//           ФУНКЦИЯ ВЕРИФИКАЦИИ
// ========================================================================
void verify()
{
    int rank_size;
    MPI_Comm_rank(main_comm, &rank_world);
    MPI_Comm_size(main_comm, &rank_size);

    int start_i =  rank_world    * (N / rank_size);
    int end_i   = (rank_world+1) * (N / rank_size);
    if (rank_world == rank_size - 1) {
        end_i = N;
    }

    double local_sum = 0.0;
    for (int i = start_i; i < end_i; i++) {
        for (int j = 0; j < N; j++) {
            local_sum += A[i][j] * (i+1) * (j+1) / (N*N);
        }
    }

    double global_sum = 0.0;
    MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, main_comm);

    if (rank_world == 0) {
        printf("S = %f\n", global_sum);
    }
}

// ========================================================================
//          ОБРАБОТЧИК ОШИБОК ДЛЯ СЦЕНАРИЯ (a)
// ========================================================================
void error_handler(MPI_Comm *pcomm, int *error_code, ...)
{
    if (has_died) {
        return;
    }
    has_died = 1;

    int wrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    fprintf(stderr, "Process %d enters error handler (code=%d)...\n", wrank, *error_code);

    // 1) Revoke: «обнуляем» текущий коммуникатор
    MPIX_Comm_revoke(*pcomm);

    // 2) Shrink: убираем упавшие процессы
    //    (после этого main_comm — новый, «усечённый»)
    MPIX_Comm_shrink(*pcomm, &main_comm);

    // Привязать обработчик ошибок и к новому коммуникатору:
    MPI_Comm_set_errhandler(main_comm, errh);

    MPI_Comm_rank(main_comm, &rank_world);
    MPI_Comm_size(main_comm, &size_world);

    // 3) Загружаем состояние из чекпоинта (или инициализируем)
    int iteration_load = 1;
    if (!load_checkpoint(&iteration_load)) {
        iteration_load = 1;
        init(/*load_from_checkpoint=*/0, &iteration_load);
    }

    it = iteration_load;

    if (rank_world == 0) {
        fprintf(stderr, "Recovered from failure. new_rank=%d, new_size=%d. "
                        "Continue from iteration %d...\n",
                rank_world, size_world, it);
    }
    // Возвращаемся в main, цикл продолжится с it+1
}

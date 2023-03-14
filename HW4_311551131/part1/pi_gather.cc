#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>

int main(int argc, char **argv)
{
    // --- DON'T TOUCH ---
    MPI_Init(&argc, &argv);
    double start_time = MPI_Wtime();
    double pi_result;
    long long int tosses = atoi(argv[1]);
    int world_rank, world_size;
    // ---
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    srand((int)time(0) * world_rank);
    long long int my_cnt = tosses / world_size, total = 0;
    if (world_rank == world_size - 1)
        my_cnt += tosses % world_size;

    for (int i = 0; i < my_cnt; i++) {
        double x = (double) rand() / (RAND_MAX + 1.0);
        double y = (double) rand() / (RAND_MAX + 1.0);
        if ((x * x + y * y) <= 1.0)
            total++;
    }
    // TODO: MPI init
    long long int cnts[world_size];
    
    // TODO: use MPI_Gather
    MPI_Gather(&total, 1, MPI_LONG, cnts, 1, MPI_LONG, 0, MPI_COMM_WORLD);

    if (world_rank == 0)
    {
        // TODO: PI result
        long long int all_total = 0;
        for (int i = 0; i < world_size; i++) {
            all_total += cnts[i];
        }
        // --- DON'T TOUCH ---
        pi_result = (double)all_total / tosses * 4;
        double end_time = MPI_Wtime();
        printf("%lf\n", pi_result);
        printf("MPI running time: %lf Seconds\n", end_time - start_time);
        // ---
    }
    
    MPI_Finalize();
    return 0;
}

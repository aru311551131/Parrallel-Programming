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
    
    // TODO: init MPI
    int stage_cnt = 1;
    long long int all_total = total;
    while(stage_cnt < world_size)
    {
        stage_cnt *= 2;
        if (world_rank % stage_cnt)
        {
            // TODO: handle workers
            MPI_Send(&all_total, 1, MPI_LONG, (world_rank / stage_cnt) * stage_cnt, 0, MPI_COMM_WORLD);
            break;
        }
        else
        {
            
            // TODO: master
            long long int tmp;          
            MPI_Status status;
            MPI_Recv(&tmp, 1, MPI_LONG, world_rank + stage_cnt / 2, 0, MPI_COMM_WORLD, &status);
            all_total += tmp; 

        }
    }
    

    if (world_rank == 0)
    {
        // TODO: process PI result
        pi_result = (double)all_total / tosses * 4;
        // --- DON'T TOUCH ---
        double end_time = MPI_Wtime();
        printf("%lf\n", pi_result);
        printf("MPI running time: %lf Seconds\n", end_time - start_time);
        // ---
    }

    MPI_Finalize();
    return 0;
}

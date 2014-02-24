#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

double get_seconds();

main(int argc, char *argv[])
{
    // variables in use
    int iproc, nproc,i;
    char host[255], message[55];
    MPI_Status status;
    char one_byte_array[1];
    one_byte_array[0] = 'a';
    double start_time;
    double average_latency;
    int number_of_experiments = 4;
    int num_procs_involved[4] = {3, 6, 9, 15};
    int proc_to_broadcast_from[4] = {5, 5, 5, 5};
    int proc_to_reduce_to[4] = {5, 5, 5, 5};
    int number_of_iterations = 5;
    int sizes_of_data[5] = {1, 512, 4096, 32768, 65536};
    double** data = malloc(sizeof(double*)*number_of_iterations);
    double** receive_buffer = malloc(sizeof(double*)*number_of_iterations);
    for(i = 0; i < number_of_iterations; i++)
    {
        data[i] = malloc(sizeof(double)*sizes_of_data[i]);
        receive_buffer[i] = malloc(sizeof(double)*sizes_of_data[i]);
    }
    
    if(iproc == 0)
    {
        for(i = 0; i < sizes_of_data[0]; i++)
        {
            (*data)[i] = i + 1;
        }
    }

    // initialization methods
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    
    // sanity check
    gethostname(host,253);
    printf("I am proc %d of %d running on %s\n", iproc, nproc,host);

    // does ping pong between procs 0 and 1 to obtain the average latency
    MPI_Barrier(MPI_COMM_WORLD);
    if(iproc == 0)
    {
        printf("%d starting ping pong (latency test)\n", iproc);
        start_time = get_seconds();
        for(i = 0; i < 1000; i++)
        {
            MPI_Send(one_byte_array, 1, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
            MPI_Recv(one_byte_array, 1, MPI_CHAR, 1, 0, MPI_COMM_WORLD, &status);
        }
        average_latency = (get_seconds() - start_time)/2000.0;
        printf("%d: average_latency = %lf\n", iproc, average_latency);
    }
    else if(iproc == 1)
    {
        for(i = 0; i < 1000; i++)
        {
            MPI_Recv(one_byte_array, 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
            MPI_Send(one_byte_array, 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
        }
    }
    
    // computes the run times of a linear broadcast
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = get_seconds();
    MPI_Bcast(data[0], sizes_of_data[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    // computes the run times of linear reduction
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = get_seconds();
    MPI_Reduce(data[0], receive_buffer[0], sizes_of_data[0], MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if(iproc == 0)
    {
        for(i = 0; i < sizes_of_data[0]; i++)
        {
            printf("%d: data[%d] = %lf", iproc, i, (*data)[i]);
            printf("%d: reduced[%d] = %lf", iproc, i, (*receive_buffer)[i]);
        }
    }
    
    

    
    
    
    
    // free allocated memory
    for(i = 0; i < number_of_iterations; i++)
    {
        free(data[i]);
        free(receive_buffer[i]);
    }
    free(data);
    free(receive_buffer);

    MPI_Finalize();
}


double get_seconds()
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ((double) tp.tv_sec + (double) tp.tv_usec * 1e-6);
}



#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

double get_seconds();
int get_dimension(int num);
void MPI_Broadcast_Hypercube(void* buffer, int count, MPI_Datatype datatype, int iproc, int nproc, int from_proc, MPI_Comm comm);
void MPI_Reduce_Hypercube(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int iproc, int nproc, int from_proc, MPI_Comm comm);

main(int argc, char *argv[])
{
    // variables in use
    int iproc, nproc,i, iteration;
    char host[255], message[55];
    MPI_Status status;
    char one_byte_array[1];
    one_byte_array[0] = 'a';
    double start_time;
    double average_latency;
    int number_of_iterations = 2;
    int sizes_of_data[5] = {1, 32768};//512};//, 4096, 32768, 65536};
    int broadcast_proc = 5; /* change this and the next to control who we reduce and broadcast to */
    int reduce_proc = 5;
    double** data = malloc(sizeof(double*)*number_of_iterations);
    double** receive_buffer = malloc(sizeof(double*)*number_of_iterations);
    for(i = 0; i < number_of_iterations; i++)
    {
        data[i] = malloc(sizeof(double)*sizes_of_data[i]);
        receive_buffer[i] = malloc(sizeof(double)*sizes_of_data[i]);
    }
    

    // initialization methods
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    
    // no point in running this if there's only one process
    if(nproc == 1)
    {
        exit(0);
    }
    
    broadcast_proc = broadcast_proc % nproc;
    // only the broadcast process will initialize the array as a sanity check
    if(iproc == broadcast_proc)
    {
        for(iteration = 0; iteration < number_of_iterations; iteration++)
        {
            for(i = 0; i < sizes_of_data[iteration]; i++)
            {
                (*(data+iteration))[i] = i + 1;
            }
        }
    }
    
    // print which machine this proc is running on
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
        printf("%d, %lf, average_latency\n", iproc, average_latency);
    }
    else if(iproc == 1)
    {
        for(i = 0; i < 1000; i++)
        {
            MPI_Recv(one_byte_array, 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
            MPI_Send(one_byte_array, 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
        }
    }
    
    
    // for different array sizes do the following broadcasts and reductions
    for(iteration = 0; iteration < number_of_iterations; iteration++)
    {
        double* buffer = *(data + iteration);
        double* recv_buf = *(receive_buffer + iteration);
        int count = sizes_of_data[iteration];
        
        // computes the run times of a linear broadcast
        MPI_Barrier(MPI_COMM_WORLD);
        start_time = get_seconds();
        MPI_Bcast(buffer, count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if(iproc == broadcast_proc) 
            printf("%d, %lf,  mpi broadcast time\n", iproc, (get_seconds()-start_time));
        
        // computes the run times of linear reduction
        MPI_Barrier(MPI_COMM_WORLD);
        start_time = get_seconds();
        MPI_Reduce(buffer, recv_buf, count, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if(iproc == reduce_proc) 
            printf("%d, %lf, mpi reduce time\n", iproc, (get_seconds()-start_time));
        
        // check that the reduction is correct
        if(iproc == reduce_proc)
        {
            for(i = 0; i < sizes_of_data[0]; i++)
            {
                if((*receive_buffer)[i] != (((*data)[i]) * nproc))
                    printf("%d: warning! reduced[%d] = %lf is incorrect", iproc, i, (*receive_buffer)[i]);
            }
        }
        
        // computes the run times of a log broadcast
        MPI_Barrier(MPI_COMM_WORLD);
        start_time = get_seconds();
        MPI_Broadcast_Hypercube(buffer, count, MPI_DOUBLE, iproc, nproc, 0, MPI_COMM_WORLD);
        if(iproc == broadcast_proc) 
            printf("%d, %lf, log broadcast time\n", iproc, (get_seconds()-start_time));
        
        // computes the run times of a log reduction
        MPI_Barrier(MPI_COMM_WORLD);
        start_time = get_seconds();
        MPI_Reduce_Hypercube(buffer, recv_buf, count, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if(iproc == reduc_proc) 
            printf("%d, %lf, log reduce time\n", iproc, (get_seconds()-start_time));
        
        if(iproc == reduce_proc)
        {
            for(i = 0; i < sizes_of_data[0]; i++)
            {
                if((*receive_buffer)[i] != (((*data)[i]) * nproc))
                    printf("%d: warning! reduced[%d] = %lf is incorrect", iproc, i, (*receive_buffer)[i]);
            }
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

void MPI_Broadcast_Hypercube(void* buffer, int count, MPI_Datatype datatype, int iproc, int nproc, int from_proc, MPI_Comm comm)
{
    if(nproc == 1) return;
    MPI_Status status; 
    iproc = (iproc - from_proc) % nproc;
    
    int i;
    int dimension = get_dimension(nproc - 1);
    int my_dimension = get_dimension(iproc);
    int bit_mask = 1;
    int my_neighbor = iproc;
    
    for(i = 0; i < dimension; i++)
    {
        my_neighbor = iproc ^ bit_mask;
        if(my_neighbor < nproc)
        {
            my_neighbor = (my_neighbor + from_proc) % nproc;
            if(iproc < bit_mask)
            {
                //printf("%d: sending to %d\n", iproc, my_neighbor);
                MPI_Send(buffer, count, datatype, my_neighbor, 0, comm);
            }
            else
            {
                if(my_dimension == get_dimension(bit_mask))
                {
                    //printf("%d: receiving a message from %d\n", iproc, my_neighbor);
                    MPI_Recv(buffer, count, datatype, my_neighbor, 0, comm, &status);
                }
                else
                {
                    // printf("%d: not doing anything\n", iproc);
                }
            }
        }
        bit_mask = bit_mask << 1;
    }
}


void MPI_Reduce_Hypercube(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int iproc, int nproc, int from_proc, MPI_Comm comm)
{
    if(nproc == 1) return; // TODO really we should just copy sendbuf to recvbuf
    MPI_Status status;
    iproc = (iproc - from_proc) % nproc;
    
    int i, j;
    int dimension = get_dimension(nproc - 1);
    int my_dimension = get_dimension(iproc);
    int bit_mask = 1;
    int my_neighbor = iproc;
    
    bit_mask = bit_mask << (dimension - 1);
    
    for(i = 0; i < dimension; i++)
    {
        my_neighbor = iproc ^ bit_mask;
        if(my_neighbor < nproc)
        {
            my_neighbor = (my_neighbor + from_proc) % nproc;
            if(my_dimension > get_dimension(bit_mask))
            {
                // printf("%d am doing nothing\n", iproc);
            }
            else
            {
                if(iproc >= bit_mask)
                {
                    // printf("%d am sending to %d\n", iproc, my_neighbor);
                    MPI_Send(sendbuf, count, datatype, my_neighbor, 0, comm);
                }
                else
                {
                    // printf("%d am receiving a message from %d\n", iproc, my_neighbor);
                    MPI_Recv(recvbuf, count, datatype, my_neighbor, 0, comm, &status);
                    if(datatype == MPI_DOUBLE)
                    {
                        // perform a sum reduce
                        if(op == MPI_SUM)
                        {
                            for(j = 0; j < count; j++)
                            {
                                *(sendbuf + j) = *(sendbuf + j) + *(recvbuf + j);
                            }
                        }
                    }
                }
            }
        }
        bit_mask = bit_mask >> 1;
    }
}


double get_seconds()
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ((double) tp.tv_sec + (double) tp.tv_usec * 1e-6);
}

int get_dimension(int num)
{
    int dimension = 0;
    for(dimension = 0; num > 0;dimension++)
    {
        num = num >> 1;
    }
    return dimension;
}

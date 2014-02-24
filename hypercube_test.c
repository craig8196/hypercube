#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DEBUG
#undef DEBUG
#define DEBUG_get_dimension
#undef DEBUG_get_dimension



double get_seconds();
int get_dimension(int num);
void MPI_Broadcast_Hypercube(int iproc, int nproc, int from_proc);
void MPI_Reduce_Hypercube(int iproc, int nproc, int from_proc);

int main(int argc, char *argv[])
{
    if(get_dimension(0) == 0 && get_dimension(1) == 1 && get_dimension(2) == 2 && get_dimension(3) == 2 &&
        get_dimension(4) == 3 && get_dimension(5) == 3 && get_dimension(6) == 3 && get_dimension(7) == 3 &&
        get_dimension(8) == 4)
    {
        printf("Passed dimension test.\n");
    }
    printf("\n");
    MPI_Broadcast_Hypercube(0, 2, 0);
    MPI_Broadcast_Hypercube(1, 2, 0);
    printf("\n");
    MPI_Broadcast_Hypercube(0, 3, 0);
    MPI_Broadcast_Hypercube(1, 3, 0);
    MPI_Broadcast_Hypercube(2, 3, 0);
    printf("\n");
    MPI_Broadcast_Hypercube(0, 7, 0);
    MPI_Broadcast_Hypercube(1, 7, 0);
    MPI_Broadcast_Hypercube(2, 7, 0);
    MPI_Broadcast_Hypercube(3, 7, 0);
    MPI_Broadcast_Hypercube(4, 7, 0);
    MPI_Broadcast_Hypercube(5, 7, 0);
    MPI_Broadcast_Hypercube(6, 7, 0);
    
    printf("\n");
    MPI_Reduce_Hypercube(0, 6, 0);
    MPI_Reduce_Hypercube(1, 6, 0);
    MPI_Reduce_Hypercube(2, 6, 0);
    MPI_Reduce_Hypercube(3, 6, 0);
    MPI_Reduce_Hypercube(4, 6, 0);
    MPI_Reduce_Hypercube(5, 6, 0);
    printf("\n");
    MPI_Reduce_Hypercube(7, 15, 0);
    printf("\n");
    MPI_Broadcast_Hypercube(6, 7, 6);
    printf("\n");
    MPI_Reduce_Hypercube(9, 15, 2);
    
    return 0;
}

/*
 * 0 = 0000
 * 1 = 0001
 * 2 = 0010
 * 3 = 0011
 * 4 = 0100
 * 5 = 0101
 * 6 = 0110
 * 7 = 0111
 * 8 = 1000
 */
void MPI_Broadcast_Hypercube(int iproc, int nproc, int from_proc)
{
    if(nproc == 1) return;
    
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
                printf("%d am sending to %d\n", iproc, my_neighbor);
            }
            else
            {
                if(my_dimension == get_dimension(bit_mask))
                {
                    printf("%d am receiving a message from %d\n", iproc, my_neighbor);
                }
                else
                {
                    printf("%d am not doing anything\n", iproc);
                }
            }
        }
        bit_mask = bit_mask << 1;
    }
}

/*
 * 0 = 0000
 * 1 = 0001
 * 2 = 0010
 * 3 = 0011
 * 4 = 0100
 * 5 = 0101
 * 6 = 0110
 * 7 = 0111
 * 8 = 1000
 */
void MPI_Reduce_Hypercube(int iproc, int nproc, int from_proc)
{
    if(nproc == 1) return;
    
    iproc = (iproc - from_proc) % nproc;
    
    int i;
    int dimension = get_dimension(nproc - 1);
    int my_dimension = get_dimension(iproc);
    int bit_mask = 1;
    int my_neighbor = iproc;
    
    bit_mask = bit_mask << (dimension - 1);
    
    for(i = 0; i < dimension; i++)
    {
        my_neighbor = iproc ^ bit_mask; //0010 0001 0011
        if(my_neighbor < nproc)
        {
            my_neighbor = (my_neighbor + from_proc) % nproc;
            if(my_dimension > get_dimension(bit_mask))
            {
                printf("%d am doing nothing\n", iproc);
            }
            else
            {
                if(iproc >= bit_mask)
                {
                    printf("%d am sending to %d\n", iproc, my_neighbor);
                }
                else
                {
                    printf("%d am receiving a message from %d\n", iproc, my_neighbor);
                }
            }
        }
        bit_mask = bit_mask >> 1;
    }
}

// num must be positive or zero
int get_dimension(int num)
{
    int dimension = 0;
#ifdef DEBUG_get_dimension
printf("Start, num = %d\n", num);
#endif
    for(dimension = 0; num > 0;dimension++)
    {
        num = num >> 1;
#ifdef DEBUG_get_dimension
printf("After a shift, num = %d\n", num);
#endif
    }
#ifdef DEBUG_get_dimension
printf("Resulting dimension = %d\n", dimension);
#endif
    return dimension;
}

double get_seconds()
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ((double) tp.tv_sec + (double) tp.tv_usec * 1e-6);
}



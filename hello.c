#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

main(int argc, char *argv[])
{
    int iproc, nproc,i;
    char host[255], message[55];
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);

    gethostname(host,253);
    printf("I am proc %d of %d running on %s\n", iproc, nproc,host);

    MPI_Finalize();
}


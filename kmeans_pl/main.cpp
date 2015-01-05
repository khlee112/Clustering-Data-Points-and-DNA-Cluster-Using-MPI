#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <getopt.h>
#include "mpi_kmeans.h"

using namespace std;

void usage(char *name) {
    printf("Usage: %s [-k <k>] [-n <numlines>] [-f <filename>] [-p|-d]\n", name);
}

int main(int argc, char *argv[]) {

    // timer
    float clocks_per_sec = 1.0 * CLOCKS_PER_SEC;    
    time_t start_p = clock();

    // MPI vars
    int numprocs, rank, namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(processor_name, &namelen);

    printf("Process %d on %s out of %d\n", rank, processor_name, numprocs);

    int linelen = 100; // hard-coded
    int whichprocess = 0, k = 10, numlines = 900;
    string filenamestring = "input/data/points.csv";

    // getopt vars
    char opt;
    while ((opt = getopt(argc, argv, "k:n:f:pdh")) != EOF) {
        switch(opt) {
            case 'k':
                k = atoi(optarg); 
                break;
            case 'n':
                numlines = atoi(optarg);
                break;
            case 'f':
                filenamestring = optarg;
                break;
            case 'p':
                whichprocess = 0;
                break;
            case 'd':
                whichprocess = 1;
                break;
            case 'h':
            default:
                usage(argv[0]);
                MPI_Finalize();
                return 0;
        }
    }

    char filename[100];
    strcpy(filename,filenamestring.c_str());

    fileinfo info = {
        .linelen = linelen,
        .numlines = numlines,
        .filename = filename
    };

    if (whichprocess==0) point_kmeans(rank, numprocs, k, info);
    else dna_kmeans(rank, numprocs, k, info);
    MPI_Finalize();

    // timer end
    time_t end_p = clock();
    float runtime_p = (end_p - start_p) / clocks_per_sec;
    printf("Rank %d Runtime: %f sec\n", rank, runtime_p);

    return 0;    
}


/**
* Utilities
*/

// displacements used in scatterv
void get_displs(int *displs, int *sendcnts, int numprocs) {
    displs[0] = 0;
    for (int i = 1; i < numprocs; i++)
        displs[i] = displs[i-1] + sendcnts[i-1];
}

// number of data to send per proc
void get_sendcnts(int *sendcnts, int numlines, int numprocs) {
    int sendcnt = numlines / numprocs;
    int leftovers = numlines % numprocs;
    for (int i = 0; i < leftovers; i++)
        sendcnts[i] = sendcnt + 1;
    for (int i = leftovers; i < numprocs; i++)
        sendcnts[i] = sendcnt;
}

// check if int n is in array a of size t
int checkdup(int n, int *a, int t) {
    for (int i = 0; i < n; i++)
        if (a[i] == t) return 1;
    return 0;
}
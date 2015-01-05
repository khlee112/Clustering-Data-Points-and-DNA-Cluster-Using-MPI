#ifndef __KMEANS_H__
#define __KMEANS_H__

typedef struct fileinfo_struct {
    const int linelen;
    const int numlines;
    const char *filename;
} fileinfo;

// Util Function
int checkdup(int n, int *a, int t);
void get_sendcnts(int *sendcnts, int numlines, int numprocs);
void get_displs(int *displs, int *sendcnts, int numprocs);

// POINT_KMEANS
int point_kmeans(int rank, int numprocs, int k, fileinfo info);

// DNA_KMEANS
int dna_kmeans(int rank, int numprocs, int k, fileinfo info);


#endif
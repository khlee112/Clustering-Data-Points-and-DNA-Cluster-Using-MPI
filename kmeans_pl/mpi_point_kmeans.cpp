#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <math.h>
#include <mpi.h>
#include "mpi_kmeans.h"

using namespace std;

typedef float point[2];

int read_data(point *points, fileinfo info);
void print_means(int k, point *means, int nonumbers);
void init_point_means(int k, point *means, int total, point *all);
void scatter_point_data(point *all, int sum, point *pts, int cnt, int np, int rank);
void addpt(point sum, point p);
float distance(point a, point b);

// mean function that initializes data and performs kmeans
int point_kmeans(int rank, int numprocs, int k, fileinfo info) {

    point all_points[info.numlines];
    point means[k];

    // allocate memory for slave use
    int remainder = info.numlines % numprocs;
    int recvcnt = info.numlines / numprocs + (rank < remainder ? 1:0);
    point points[recvcnt];

    // read data and initialize
    if (!rank) { // if master (rank 0)
        read_data(all_points, info);
        init_point_means(k, means, info.numlines, all_points);
    }
    MPI_Bcast(means, 2 * k, MPI_FLOAT, 0, MPI_COMM_WORLD);

    // scatter the data
    scatter_point_data(all_points, info.numlines, points, recvcnt, numprocs, rank);

    // iterate
    int changed;
    do {
        changed = 0;
        point partial_sum[k];
        int counts[k];
        memset(partial_sum, 0, sizeof(partial_sum));
        memset(counts, 0, sizeof(counts));

        // determine assignments and sum partial means
        for (int i = 0; i < recvcnt; i++) {
            int closest = 0;
            float dist = distance(means[0], points[i]);
            float distj;
            for (int j = 1; j < k; j++) {
                distj = distance(means[j], points[i]);
                if (distj < dist) {
                    closest = j;
                    dist = distj;
                }
            }
            addpt(partial_sum[closest], points[i]);
            counts[closest]++;
        }
        
        // reduce clusters to centroids
        int clusters[k];
        point newmeans[k];

        MPI_Allreduce(counts, clusters, k, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(partial_sum, newmeans, 2*k, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

        for (int i = 0; i < k; i++) {
            if (!clusters[i]) continue;

            point mean;
            memcpy(mean, means[i], sizeof(mean));

            float x = newmeans[i][0] / clusters[i];
            float y = newmeans[i][1] / clusters[i];

            if (x != mean[0] || y != mean[1]) {
                changed = 1;
                mean[0] = x;
                mean[1] = y;
                memcpy(means[i], mean, sizeof(mean));
            }
        }
        // print mean
        //if (!rank) print_means(k, means, 1);

    } while (changed);

    // print mean
    if (!rank) print_means(k, means, 1);
}

void addpt(point sum, point p) {
    sum[0] += p[0];
    sum[1] += p[1];
}

void print_means(int k, point *means, int nonumbers) {
    printf("updated means:\n");
    for (int i = 0; i < k; i++) {
        if (!nonumbers) printf("%d: ", i);
        printf("(%f, %f)\n", means[i][0], means[i][1]);
    }
    printf("\n");
}

// initialize random means by choosing points from data
void init_point_means(int k, point *means, int total, point *all) {
    srand(1); // arbitrary but non-random seed for data analysis
    int meanis[k];
    memset(meanis, -1, k);

    for (int i = 0; i < k; i++) {
        int index = rand() % total;
        if (k < total)
            while (checkdup(k, meanis, index))
                index = rand() % total;
            
        meanis[i] = index;
        memcpy(means[i], all[index], sizeof(means[i]));
    }
}

int read_data(point *points, fileinfo info) {
    FILE *fp;
    if (!(fp = fopen(info.filename, "r")))
        return -1;

    char *line = NULL;
    size_t len = 0; 
    size_t read;
    int numpoints = 0;
    for (int lineno = 0; lineno < info.numlines; lineno++) {
        read = getline(&line, &len, fp);
        point p = {
            atof(line),
            atof(strchr(line, ',') + 1)
        };
        memcpy(points[numpoints++], p, sizeof(p));
    }
    if (line) free(line);
    fclose(fp);
    return 0;
}

// each proc gets a fraction of the total data points
void scatter_point_data(point *all, int sum, point *pts, int cnt, int np, int rank) {
    int sendcnts[np], displs[np];
    get_sendcnts(sendcnts, sum, np);
    for (int i = 0; i < np; i++)
        sendcnts[i] *= 2;
    get_displs(displs, sendcnts, np);
    MPI_Scatterv(all, sendcnts, displs, MPI_FLOAT, pts, 2*cnt, MPI_FLOAT,
            0, MPI_COMM_WORLD);
}

// euclidean distance between two points
float distance(point a, point b) {
    return sqrt(pow(a[0] - b[0], 2) + pow(a[1] - b[1], 2));
}

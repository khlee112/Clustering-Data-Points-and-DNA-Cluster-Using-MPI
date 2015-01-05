#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/stat.h>
#include <math.h>
#include <mpi.h>
#include <assert.h>
#include "mpi_kmeans.h"

using namespace std;

float distance(char* a, char* b);
void normalize(int *arr, int numprocs);
void scatter_dna_data(char *all, int total, char *recvbuf, int revcnt, int np, int rank);
int read_data(char *recvbuf, fileinfo info);
void init_dna_means(int k, char* means, int total, char *all);

int strand_len;

int dna_kmeans(int rank, int numprocs, int k, fileinfo info) {

    strand_len = info.linelen + 1; // end with '\0'
    
    // allocate memory for master use
    char* all_strands = (char*)malloc(sizeof(char) * info.numlines * strand_len);
    char* means = (char*)malloc(sizeof(char) * k * strand_len);
    memset(all_strands, 0, sizeof(all_strands));
    memset(means, 0, sizeof(means));

    // allocate memory for slave use
    int remainder = info.numlines % numprocs;
    int recvcnt = info.numlines / numprocs + (rank < remainder ? 1:0);
    char* recvbuf = (char*)malloc(sizeof(char) * recvcnt*strand_len);

    // read data and initialize
    if (!rank) { // if master (rank 0)
        read_data(all_strands, info);
        init_dna_means(k, means, info.numlines, all_strands);
    }
    MPI_Bcast(means, k*strand_len, MPI_CHAR, 0, MPI_COMM_WORLD);

    // scatter data
    scatter_dna_data(all_strands, info.numlines, recvbuf, recvcnt, numprocs, rank);

    // iterate
    int changed;
    do {
        changed = 0;
        char *kmeans[k][recvcnt];
        int counts[k];
        memset(kmeans, 0, sizeof(kmeans));
        memset(counts, 0, sizeof(counts));

        // calculate recvcnt distances
        for (int i = 0; i < recvcnt; i++) {
            int closest = 0;
            float dist = distance(means, &recvbuf[i*strand_len]);
            float jdist;
            for (int j = 1; j < k; j++) {
                jdist = distance(&means[j*strand_len], &recvbuf[i*strand_len]);
                if (jdist < dist) {
                    closest = j;
                    dist = jdist;
                }
            }
            kmeans[closest][counts[closest]++] = recvbuf + i*strand_len;
        }
        
        // reduce clusters to centroids
        for (int i = 0; i < k; i++) {
            char newmean[strand_len];
            memset(newmean, 0, sizeof(newmean));
            // for each position
            for (int t = 0; t < strand_len - 1; t++) {
                int freq[4];
                memset(freq, 0, sizeof(freq));

                for (int j = 0; j < counts[i]; j++) {
                    switch(kmeans[i][j][t]) {
                        case 'A':
                            freq[0]++;
                            break;
                        case 'C':
                            freq[1]++;
                            break;
                        case 'G':
                            freq[2]++;
                            break;
                        case 'T':
                            freq[3]++;
                            break;
                        default:
                            assert(0);
                    }
                }

                int countA, countC, countG, countT;
                MPI_Reduce(freq, &countA, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
                MPI_Reduce(&freq[1], &countC, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
                MPI_Reduce(&freq[2], &countG, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
                MPI_Reduce(&freq[3], &countT, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

                if (!rank) {
                    newmean[t] = 'A';
                    int cnt = countA;
                    if (countC > cnt) {
                        newmean[t] = 'C';
                        cnt = countC;
                    }
                    if (countG > cnt) {
                        newmean[t] = 'G';
                        cnt = countG;
                    }
                    if (countT > cnt) {
                        newmean[t] = 'T';
                        cnt = countT;
                    }
                }
            }

            if (!rank) {
                if (strncmp(means+i*strand_len, newmean, strand_len-1)) {
                    changed = 1;
                    // Comment Out
                    /*
                    printf("updating mean[%d]: %s-->%s\n",
                        i, means+i*strand_len, newmean);
                    */
                    memcpy(means+i*strand_len, newmean, strand_len-1);
                }
            }
        }

        MPI_Bcast(&changed, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(means, k*strand_len, MPI_CHAR, 0, MPI_COMM_WORLD);

        if (!rank) {
            /*
            // Comment Out
            for (int i = 0; i < k; i++)
                printf("update mean %d: %s\n", i, means+i*strand_len);
            */
        }

    } while (changed);

    if (!rank) printf("process done.\n");
    free(all_strands);
    free(means);
    free(recvbuf);
}

void init_dna_means(int k, char* means, int total, char *all) {
    srand(time(NULL));
    int meanis[k];
    memset(meanis, -1, k);

    for (int i = 0; i < k; i++) {
        int index = rand() % total;
        if (k < total) {
            while (checkdup(k, meanis, index))
                index = rand() % total;
        }    
        meanis[i] = index;
        memcpy(means+i*strand_len, all+index*strand_len, 
            (strand_len)*sizeof(char));
    }
}

int read_data(char *recvbuf, fileinfo info) {

    FILE *fp;
    if (!(fp = fopen(info.filename, "r"))) return -1;

    char line[info.linelen+1];
    int num = 0;
    while (fgets(line, sizeof(line), fp)!= NULL) {
       if (line[strlen(line)] == '\n') line[strlen(line)] = '\0';
       if (strlen(line) < info.linelen) continue;
       strncpy(recvbuf + num, line, strlen(line));
       num += strlen(line);
       recvbuf[num++] = '\0';
    }
    fclose(fp);
    return 0;
}

void scatter_dna_data(char *all, int total, char *recvbuf, 
    int revcnt, int np, int rank) {
    
    int sendcnts[np], displs[np];
    get_sendcnts(sendcnts, total, np);
    get_displs(displs, sendcnts, np);
    // normalize the size
    normalize(sendcnts, np);
    normalize(displs, np);
    revcnt *= strand_len;
    MPI_Scatterv(all, sendcnts, displs, MPI_CHAR, recvbuf, revcnt,
                MPI_CHAR, 0, MPI_COMM_WORLD);
}

void normalize(int *arr, int numprocs) {
    for (int i=0; i < numprocs; i++) arr[i] *= strand_len;
}

float distance(char* a, char* b) {
    float diff = 0;
    for (int i = 0; i < strand_len; i++)
        if (a[i] != b[i]) diff++;
    return diff;
}


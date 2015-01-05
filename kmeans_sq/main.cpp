#include "dna_strand.h"
#include "point.h"
#include "kmeans_seq.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string>
#include <cstring>
#include <random>

using namespace std;

void usage(char *name) {
    printf("Usage: %s [-k <k>][-f <filename>] [-p|-d]\n", name);
}

int main(int argc, char *argv[]) {
	
	// timer
	float clocks_per_sec = 1.0 * CLOCKS_PER_SEC;
	time_t start_p = clock();
	
    int whichprocess = 0, k = 10;
    string filenamestring = "input/data/points.csv";

	// getopt vars
	char opt;
    while ((opt = getopt(argc, argv, "k:n:f:pdh")) != EOF) {
        switch(opt) {
            case 'k':
                k = atoi(optarg); 
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
                return 0;
        }
    }

    char filename[100];
    strcpy(filename,filenamestring.c_str());

	int tolerance = 0; // allowed label changes for convergence

	if (whichprocess==0) {
		Kmeans_Seq <Point> kmeans_on_point(filename, k, tolerance);
    	kmeans_on_point.kmeans();
		kmeans_on_point.print_centroids();
	} else {
		Kmeans_Seq <DNA_Strand> kmeans_on_dna_strands(filename, k, tolerance);
    	kmeans_on_dna_strands.kmeans();
    	kmeans_on_dna_strands.print_centroids();
	}

	// timer end
    time_t end_p = clock();
    float runtime_p = (end_p - start_p) / clocks_per_sec;
    printf("Runtime: %f sec\n", runtime_p);
	return 0;
}

/**
* Utilities
*/

// check if int n is in array a of size t
int checkdup(int n, int *a, int t) {
    for (int i = 0; i < n; i++)
        if (a[i] == t) return 1;
    return 0;
}

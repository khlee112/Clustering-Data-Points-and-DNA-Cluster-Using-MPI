#ifndef _KMEANS_SEQ_H
#define _KMEANS_SEQ_H

#include <vector>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <iostream>
#include <sys/stat.h>
#include <fstream>
#include <math.h>

using namespace std;

int checkdup(int n, int *a, int t);

template<class Data_Type>
class Kmeans_Seq {
public:
	vector<int> label; // cluster label on data points
	vector<Data_Type> centroids; // centroids of clusters

	Kmeans_Seq(string, int, int);
	void read_data(string filename);
	void kmeans();
	void print_centroids();

private:
	int N; // N data
	int K; // K cluster
	int tolerance; // converange threshold
	vector<Data_Type> data;

	void initialize();
};

template <class Data_Type> Kmeans_Seq<Data_Type>::Kmeans_Seq(string filename,
    int _K, int _tolerance) {
	K = _K;
	tolerance = _tolerance;
	read_data(filename);
	label.resize(N);
	centroids.resize(K);
}

template <class Data_Type> void Kmeans_Seq<Data_Type>::read_data(string filename) {
    string line;
    Data_Type d;
	ifstream infile(filename.c_str());
	if (infile.is_open()) {
	    while (getline(infile, line)) {
	    	//cout << line << "\n";
	        d = Data_Type(line);
	        data.push_back(d);
	    }
	    infile.close();
	}
    N = data.size();
}

/**
* assign initial cluster centroids
*/
template <class Data_Type> void Kmeans_Seq<Data_Type>::initialize() {
    srand(1); // arbitrary but non-random seed for data analysis
	int meanis[K];
	memset(meanis, -1, K);

	fill(label.begin(), label.end(), 0);
	for (int i=0; i<K; i++) {
		int index = rand() % N;
		if (K < N)
            while (checkdup(K, meanis, index))
                index = rand() % N;

		centroids[i] = data[index]; // pick first K data points
	}
}

/**
* compute new centroids
*/
template <class Data_Type> void Kmeans_Seq<Data_Type>::kmeans() {
	initialize();
	int label_change_count;

	do {
		vector<vector<Data_Type> > clusters (K); // K x (number of points in cluster) vector
		label_change_count = 0;

		// Clustering
		for (int i=0; i<N; i++) {
			float min_distance = data[i].get_distance(centroids[0]);
			int cluster_id = 0;
			for (int j=1; j<K; j++) {
				int distance = data[i].get_distance(centroids[j]);
				if (distance < min_distance) {
					cluster_id = j;
					min_distance = distance;
				}
			}
			if (label[i] != cluster_id) {
				label_change_count++;
				label[i] = cluster_id;
			}
			clusters[cluster_id].push_back(data[i]);
		}
		
		// Compute mean point of each class
		for (int j=0; j<K; j++) {
			if (clusters[j].size() > 0) centroids[j] = Data_Type::mean(clusters[j]);
		}
		// cout << "label_change_count: " << label_change_count << "\n";

	} while (label_change_count > tolerance);
}

template <class Data_Type> void Kmeans_Seq<Data_Type>::print_centroids() {
	for (int i=0; i<K; i++) {
		cout<<"("<<centroids[i].to_string()<<")\n";
	}
}

#endif

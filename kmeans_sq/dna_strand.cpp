#include "kmeans_seq.h"
#include "dna_strand.h"
#include <math.h>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <map>

using namespace std;

/**
* Implement DNA_Strand class defined in dna_strand.h
*/

DNA_Strand::DNA_Strand(string s) {
	strand = s;
}

float DNA_Strand::get_distance(DNA_Strand s) {
	float distance = 0.0;
	for (int i=0; i<s.strand.size(); i++) {
		if (s.strand[i] != strand[i]) distance += 1;
	}
	return distance;
}

string DNA_Strand::to_string() {
	return strand; 
}

/**
* compute average DNA strand
*/
DNA_Strand DNA_Strand::mean(vector<DNA_Strand> strands) {
	
	DNA_Strand mean_strand;
	mean_strand.strand = "";
	int N = strands.size();
	int length = strands[0].strand.length();
	vector<map<char, float> > buckets(length);

	for (int i=0; i<length; i++) {
		buckets[i]['A'] = 0;
		buckets[i]['C'] = 0;
		buckets[i]['G'] = 0;
		buckets[i]['T'] = 0;
	}
	
	for (int j = 0; j < N; j++) {
		for (int i=0; i < length; i++) {
			buckets[i][strands[j].strand[i]] += 1;
		}
	}

	for (int i=0; i<length; i++) {

		float max = buckets[i]['A'];
		char key = 'A';
		
		if (buckets[i]['C'] > max) {
			max = buckets[i]['C'];
			key = 'C';
		} 
		
		if (buckets[i]['G'] > max) {
			max = buckets[i]['G'];
			key = 'G';
		}
		
		if (buckets[i]['T'] > max) {
			max = buckets[i]['T'];
			key = 'T';
		}

		mean_strand.strand.append(1, key);
	}

	return mean_strand;
}
#ifndef _DNA_STRAND_H
#define _DNA_STRAND_H

#include <vector>
#include <string>

using namespace std;

/**
* Represent DNA strand composed of 'A''C''G''T
*/
class DNA_Strand {
public:
	string strand;
	DNA_Strand() {};
	DNA_Strand(string);
	float get_distance(DNA_Strand);
	string to_string();
	static DNA_Strand mean(vector<DNA_Strand>);
};

#endif
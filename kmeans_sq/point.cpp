#include "kmeans_seq.h"
#include "point.h"
#include <math.h>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>

using namespace std;

/**
* Implement Point class defined in point.h
*/

Point::Point(string s) {
	x = atof(s.c_str());
	y = atof(strchr(s.c_str(), ',') + 1);
}

Point::Point(float _x, float _y) {
	x = _x;
	y = _y;
}

float Point::get_distance(Point p) {
	return sqrt((x-p.x)*(x-p.x)+(y-p.y)*(y-p.y));
}

string Point::to_string() {
	ostringstream stream;
  	stream << x << ", " << y;
	return stream.str(); 
}

Point Point::mean(vector<Point> points) {
	int N = points.size();
	Point point;
	point.x = 0;
	point.y = 0;
	
	for (int i = 0; i < N ; i++) {
		point.x += points[i].x/N;
		point.y += points[i].y/N;
	}

	return point;
}
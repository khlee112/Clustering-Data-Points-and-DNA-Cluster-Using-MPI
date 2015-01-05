#ifndef _POINT_H
#define _POINT_H

#include <vector>
#include <string>

using namespace std;

/**
* Represent 2D data point(x, y)
*/
class Point {
public:
	float x,y;
	Point() {};
	Point(string);
	Point(float, float);
	float get_distance(Point);
	string to_string();
	static Point mean(vector<Point>);
};

#endif
#include <time.h>

void randomize(){
	//srand((time_t) time(NULL));  
}

double norm(Point v){
	return sqrt( pow(v.x(),2) + pow(v.y(),2) + pow(v.z(),2));
}

Point unitNorm(Point &v){
	double n = norm(v);
	return Point(v.x()/n, v.y()/n, v.z()/n );
}

double norm(Vector &v){
	return sqrt( v*v );
}

Vector unit_vector(Vector &v){
	double n = norm(v);
	return v/n;
}

bool IsNumber(double x) {
    return (x == x); 
}
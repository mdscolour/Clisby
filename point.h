#pragma once
#include "global.h"

class Point
{
public:
	double x, y, z;
	//everything just inline
	Point(){}
	Point(const char* key){ if (key == "zero")this->zero(); }
	Point(double a, double b, double c) :x(a), y(b), z(c){}
	//GPoint(const Sphere& p) :x(p.x), y(p.y), z(p.z){printf("aaa done once\n");}
	//GPoint(const GPoint& p) :x(p.x), y(p.y), z(p.z){}
	Point& operator=(const Point& p){ x = p.x; y = p.y; z = p.z; return(*this); }

	double coord_x(){ return(x); }
	double coord_y(){ return(y); }
	double coord_z(){ return(z); }
	Point operator+(Point p){ return Point(x + p.x, y + p.y, z + p.z); }
	Point operator-(Point p){ return Point(x - p.x, y - p.y, z - p.z); }
	Point operator*(double c){ return Point(x*c, y*c, z*c); }
	Point operator/(double c){ return Point(x/c, y/c, z/c); }
	Point& operator+=(Point p){ x += p.x; y += p.y; z += p.z; return *this; }
	Point& operator-=(Point p){ x -= p.x; y -= p.y; z -= p.z; return *this; }
	Point& operator*=(double c){ x *= c; y *= c; z *= c; return *this; }
	Point& operator/=(double c){ x /= c; y /= c; z /= c; return *this; }

	//bool operator==(GPoint p){ return((abs(p.x - x)<1e-8 && abs(p.y - y)<1e-8 && abs(p.z - z)<1e-8)); }
	//bool close(GPoint p, double preci = 1.e-5){ return((abs(p.x - x)<preci && abs(p.y - y)<preci && abs(p.z - z)<preci)); }

	void point_assign(double xx, double yy, double zz) { x = xx; y = yy; z = zz; }
	double dot(Point p){ return(x*p.x+y*p.y+z*p.z); }

	void zero()	{x = 0;y = 0;z = 0;}
	//void assign(T xx, T yy, T zz) { x = xx; y = yy; z = zz; }

	double norm(){return(double(sqrt(x*x + y*y + z*z)));}
};

// 3x3 matrix
class Matrix
{
public:
	Point row1;
	Point row2;
	Point row3;
	Matrix(){}
	Matrix(const char* key){ if (key == "identity")this->identity(); }
	Matrix(Point a, Point b, Point c) :row1(a), row2(b), row3(c){}
	Matrix(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3) :row1(x1, y1, z1), row2(x2, y2, z2), row3(x3, y3, z3){}
	//Matrix(RMatrix x) :row1(x.row1), row2(x.row2), row3(x.row3){}

	Point col1(){ return Point(row1.coord_x(), row2.coord_x(), row3.coord_x()); }
	Point col2(){ return Point(row1.coord_y(), row2.coord_y(), row3.coord_y()); }
	Point col3(){ return Point(row1.coord_z(), row2.coord_z(), row3.coord_z()); }
	Matrix& operator=(const Matrix& x){ row1 = x.row1;  row2 = x.row2; row3 = x.row3; return(*this); }
	//bool operator==(Matrix m){ return((row1 == m.row1 && row2 == m.row2 && row3 == m.row3)); }
	//bool close(Matrix m, double preci = 0.01){ return((row1.close(m.row1, preci) && row2.close(m.row2, preci) && row3.close(m.row3, preci))); }
	Matrix inv()
	{
		double determinant = row1.x*(row2.y*row3.z - row3.y*row2.z)
			- row1.y*(row2.x*row3.z - row2.z*row3.x)
			+ row1.z*(row2.x*row3.y - row2.y*row3.x);
		double invdet = 1 / determinant;
		return Matrix(
			(row2.y*row3.z - row3.y*row2.z)*invdet,
			-(row1.y*row3.z - row1.z*row3.y)*invdet,
			(row1.y*row2.z - row1.z*row2.y)*invdet,
			-(row2.x*row3.z - row2.z*row3.x)*invdet,
			(row1.x*row3.z - row1.z*row3.x)*invdet,
			-(row1.x*row2.z - row2.x*row1.z)*invdet,
			(row2.x*row3.y - row3.x*row2.y)*invdet,
			-(row1.x*row3.y - row3.x*row1.y)*invdet,
			(row1.x*row2.y - row2.x*row1.y)*invdet);
	}

	//void print(FILE *fptr){ row1.print(fptr); row2.print(fptr); row3.print(fptr); }
	Matrix& identity(){ row1 = Point(1, 0, 0); row2 = Point(0, 1, 0); row3 = Point(0, 0, 1); return(*this); }

	Point dot(Point v){ return Point(row1.dot(v), row2.dot(v), row3.dot(v)); }		
	Matrix dot(Matrix x){ return Matrix(row1.dot(x.col1()), row1.dot(x.col2()), row1.dot(x.col3()), row2.dot(x.col1()), row2.dot(x.col2()), row2.dot(x.col3()), row3.dot(x.col1()), row3.dot(x.col2()), row3.dot(x.col3())); }
	//Sphere dot(Sphere v){return GPoint(row1.dot(v),row2.dot(v),row3.dot(v));}
};

#include "point.h"

class Box
{
public:
	//double a.x,b.x,a.y,b.y,a.z,b.z;
	Point a,b;
	Box(){}
	Box(Point aa,Point bb){a=aa;b=bb;}
	Box(double n1,double n2,double n3,double n4,double n5,double n6)
	{a.x=n1,b.x=n2,a.y=n3,b.y=n4,a.z=n5,b.z=n6;}
	Box bound(Box bb){return Box(min(this->a.x,bb.a.x),
		                         max(this->b.x,bb.b.x),
								 min(this->a.y,bb.a.y),
								 max(this->b.y,bb.b.y),
								 min(this->a.z,bb.a.z),
								 max(this->b.z,bb.b.z));}
	Box intersec(Box bb){return Box(max(this->a.x,bb.a.x),
		                         min(this->b.x,bb.b.x),
		                         max(this->a.y,bb.a.y),
		                         min(this->b.y,bb.b.y),
		                         max(this->a.z,bb.a.z),
		                         min(this->b.z,bb.b.z));}
	Box& operator+(Point x){a+=x;b+=x;return(*this);}
	bool notempty(){if(a.x>b.x)return false;if(a.y>b.y)return false;if(a.z>b.z)return false; return true;}
	Box apply(Matrix m,Point p)
	{
		Point v2=m.dot(Point(a.x,a.y,b.z)),
			v3=m.dot(Point(a.x,b.y,a.z)),
			v4=m.dot(Point(b.x,a.y,a.z)),
			v5=m.dot(Point(b.x,b.y,a.z)),
			v6=m.dot(Point(a.x,b.y,b.z)),
			v7=m.dot(Point(b.x,a.y,b.z)),
			v1=m.dot(a),
			v8=m.dot(b);
		return Box(min(min(min(v1.x,v2.x),min(v3.x,v4.x)),min(min(v5.x,v6.x),min(v7.x,v8.x))),
			max(max(max(v1.x,v2.x),max(v3.x,v4.x)),max(max(v5.x,v6.x),max(v7.x,v8.x))),
			min(min(min(v1.y,v2.y),min(v3.y,v4.y)),min(min(v5.y,v6.y),min(v7.y,v8.y))),
			max(max(max(v1.y,v2.y),max(v3.y,v4.y)),max(max(v5.y,v6.y),max(v7.y,v8.y))),
			min(min(min(v1.z,v2.z),min(v3.z,v4.z)),min(min(v5.z,v6.z),min(v7.z,v8.z))),
			max(max(max(v1.z,v2.z),max(v3.z,v4.z)),max(max(v5.z,v6.z),max(v7.z,v8.z))))+p;
	}
	Box apply(Matrix m)
	{
		Point v2=m.dot(Point(a.x,a.y,b.z)),
			v3=m.dot(Point(a.x,b.y,a.z)),
			v4=m.dot(Point(b.x,a.y,a.z)),
			v5=m.dot(Point(b.x,b.y,a.z)),
			v6=m.dot(Point(a.x,b.y,b.z)),
			v7=m.dot(Point(b.x,a.y,b.z)),
			v1=m.dot(a),
			v8=m.dot(b);
		return Box(min(min(min(v1.x,v2.x),min(v3.x,v4.x)),min(min(v5.x,v6.x),min(v7.x,v8.x))),
			       max(max(max(v1.x,v2.x),max(v3.x,v4.x)),max(max(v5.x,v6.x),max(v7.x,v8.x))),
				   min(min(min(v1.y,v2.y),min(v3.y,v4.y)),min(min(v5.y,v6.y),min(v7.y,v8.y))),
				   max(max(max(v1.y,v2.y),max(v3.y,v4.y)),max(max(v5.y,v6.y),max(v7.y,v8.y))),
				   min(min(min(v1.z,v2.z),min(v3.z,v4.z)),min(min(v5.z,v6.z),min(v7.z,v8.z))),
				   max(max(max(v1.z,v2.z),max(v3.z,v4.z)),max(max(v5.z,v6.z),max(v7.z,v8.z))));
	}
};

class tree
{
public:
	tree(){}
	tree *wl,*wr,*wp;
	int n;
	Box B;
	Point Xe,X;
	double X2;
	Matrix q;
	bool full(){return (wl!=0 && wr!=0);}
};

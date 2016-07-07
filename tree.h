
// Only for small program it doesn't matter whether everything is inline

#include "global.h"
#include "point.h"

class box
{
public:
	box(){};
	//box(Sphere aa, Sphere bb){a=aa;b=bb;}

	//box apply(Matrix m){ return box(m.dot(a),m.dot(b)); }
	Sphere a, b;
};

class node
{
public:
	//bool IsLeftChild(node* lc){return lc==this;}
	//bool IsRightChild(){return father->rchild==this;}
	void print(){printf("%d %d %d %d %d \n",(int)isleaf,npoint,(int)lchild,(int)rchild,(int)father);B.a.print(stdout);B.b.print(stdout);}

	node():isleaf(false),p("identity"){}

	node *lchild, *rchild, *father; // right child is alway lchild +1
	bool isleaf;
	int npoint;
		
/*	int id;*/
	Matrix p;
	GPoint<double> Xe;  //end point
	//GPoint<int> Xs;  //start point alway (0,1)
	box B;
};

class tree
{
public:
	tree(int tsteps); //n is steps

	int twopowof(int power){int temp = 1; for(int i=0;i<power;i++) temp*=2; return temp;}
	int factorial(int n){if(n==1|| n==0)return 1;else return n*factorial(n-1);}

	void initialize(int tsteps);

	void LineInit();


	class Proposal
	{
	public:
		Matrix RandMatrix();
		int pivot_loc;
		Matrix op;
		Matrix invop;

		Proposal(int n):pivot_loc(int(n*RNG_NAME())){op = RandMatrix();invop=op.inv();}
	};

	void Record();

	void WriteDown();

	Sphere GetStepi(int i);

	int run(int discard=0, int outer=1, int inner=10000);

	int GoNSteps(int num);

	bool PivotAttempt();

	inline void DoPivot(int pivot_loc,Matrix p);

	bool CheckPivot(Proposal prop);

	bool Intersect(node* lc,node* rc,GPoint<double>& relXe,Matrix& p);

	node* Ind2Pn(int i)
	{
		int pid = zeroid + i/2;
		if(pid>nnode-1) pid -= (nsteps+1)/2;
		return nodelist+pid;
	}

	Sphere SphereMerge(Sphere& a, Sphere& b)
	{
		GPoint<double> c = b.center()-a.center();
		double r = (c.norm()+a.r+b.r)/2.0;
		Sphere rc = a+c*((r-a.r)/c.norm());
		return Sphere(rc.x, rc.y, rc.z, r, b.k);
	}

	bool SphereIntersec(Sphere& a, Sphere& b) {return ((a.center()-b.center()).norm()> (a.r+b.r))? false:true;}

	box BoxMerge(box& b1,box& b2,GPoint<double>& b1Xe,Matrix& p)
	{
		box t; 
		t.a=SphereMerge(b1.a,b1.b); 
		t.b=p.dot(SphereMerge(b2.a,b2.b))+b1Xe;
		return t;
	}
	bool BoxIntersec(box& b1,box& b2,GPoint<double>& b1Xe,Matrix& p)  // true is intersected
	{
		if(SphereIntersec(b1.b,p.dot(b2.a)+b1Xe)) return true;
		if(SphereIntersec(b1.b,p.dot(b2.b)+b1Xe)) return true;
		if(SphereIntersec(b1.a,p.dot(b2.b)+b1Xe)) return true;
		if(SphereIntersec(b1.a,p.dot(b2.a)+b1Xe)) return true;
		return false;
	}

	void UpdateNodeData(node* pn)
	{
		pn->B = BoxMerge(pn->lchild->B,pn->rchild->B,pn->lchild->Xe,pn->p);
		pn->Xe = pn->p.dot(pn->rchild->Xe) + pn->lchild->Xe;
	}
	
	int zeroid;
	int ncheck;
	const char *data_fname, *final_name;
	node* nodelist;
	//node* ptemp;
	int nsteps;
	int deep;
	int nnode;
	int leafstartid;
};
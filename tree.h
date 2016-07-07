
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

	Sphere SphereMerge(Sphere& a, Sphere& b);

	bool SphereIntersec(Sphere& a, Sphere& b);

	box BoxMerge(box& b1,box& b2,GPoint<double>& b1Xe,Matrix& p);

	bool BoxIntersec(box& b1,box& b2,GPoint<double>& b1Xe,Matrix& p);  // true is intersected

	void UpdateNodeData(node* pn);
	
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
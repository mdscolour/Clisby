
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
	tree(int tsteps):data_fname("data"),final_name("FinalWalk"),ncheck(0) //n is steps
	{
#ifdef linux
		struct timeval tpstart;
		gettimeofday(&tpstart, NULL);
		srand48(tpstart.tv_usec);
		//srand48(floor(mytime()));
#endif
#ifdef _WIN32
		srand(unsigned int(time(NULL)));
#endif

		if ( (tsteps+1) % 2 != 0 ) 
		{
			printf("odd number of site not yet supported\n");
			tsteps++;
		}

		nsteps = tsteps;

		deep = 1; // root is deep 0
		while (twopowof(deep)<nsteps+1) deep++;
		deep--; // last level node of 1 point is not count 

		nnode = 0;
		for (int i=0;i<deep;i++) nnode += twopowof(i);
		nnode += (nsteps+1)-twopowof(deep);

		leafstartid = nnode - (nsteps+1)/2;

		nodelist = new node[nnode];

		int* level = new int[nnode];
		int* indinlevel = new int[nnode];
		int tlevel = 0;
		int tindinlevel = 0;
		for (int i=0;i<nnode;i++)
		{
			//nodelist[i].id = i;
			level[i] = tlevel;
			indinlevel[i] = tindinlevel;
			tindinlevel ++;
			if(tindinlevel>twopowof(tlevel)-1) {tindinlevel = 0; tlevel++;}			
			if(i>=leafstartid) {nodelist[i].isleaf = true; nodelist[i].npoint = 2;}
		}

		nodelist[0].father = 0;
		for (int i=leafstartid-1;i>=0;i--)
		{
			int lc = i+twopowof(level[i])+(indinlevel[i]);
			nodelist[i].npoint = nodelist[lc].npoint + nodelist[lc+1].npoint;
			nodelist[i].lchild = nodelist+lc;
			nodelist[i].rchild = nodelist+lc+1;
			nodelist[lc].father = nodelist+i;
			nodelist[lc+1].father = nodelist+i;
		}

		LineInit();

		delete[] level;
		delete[] indinlevel;

// 		for (int i=0;i<nnode;i++)
// 		{
// 			nodelist[i].print();
// 		}
	}

	class Proposal
	{
	public:
		int pivot_loc;
		Matrix op;
		Matrix invop;

		Proposal(int n) :pivot_loc(int(n*RNG_NAME())){op.rand();invop=op.inv();}
	};

	void Record()
	{
		FILE *fptr;
		// record the "data"
		double endnorm = nodelist[0].Xe.norm();
		fptr = fopen(data_fname, "a");
		fprintf(fptr,"%14.10f    %14.10f \n",123.0,endnorm*endnorm);
		fclose(fptr);
	}

	void WriteDown()
	{
		FILE *fptr;
		fptr = fopen(final_name, "w");

		int i;
		fprintf(fptr, "%d  %lf  \n", nsteps, -1.0);
		for (i = 0; i <= nsteps; i++) GetStepi(i).print(fptr);

		fclose(fptr);
	}

	Sphere GetStepi(int i)
	{
		Sphere t;
		node* pn;
		if (i % 2 == 1)
		{
			pn = &nodelist[nnode - (nsteps-i)/2 -1];
			t = pn->B.b;
		}
		else
		{
			pn = &nodelist[nnode - (nsteps-i-1)/2 -1];
			t = pn->B.a;
		}
		while (pn != nodelist)
		{
			if(pn==nodelist) printf("GetStepi error\n");
			if(pn->father->rchild==pn) { t = pn->father->p.dot(t)+pn->father->lchild->Xe;}
			pn = pn->father;
		}
		return t;
	}

	int run(int discard=0, int outer=1, int inner=1000)
	{
		int accept = 0;
		accept += GoNSteps(discard);
		for (int i=0;i<outer;i++)
		{
			accept += GoNSteps(inner);
			Record();
			printf("%lf \n",ncheck/(double)inner);
			ncheck = 0;
		}
		WriteDown();
		return accept;
	}

	int GoNSteps(int num)
	{
		int accept = 0;
		for (int i=0;i<num;i++)
		{
			if(PivotAttempt()) accept++;
		}
		return accept;
	}

	void DoPivot(int pivot_loc,Matrix p)
	{
		node* pn;
		if (pivot_loc % 2 == 1)
		{
			pn = &nodelist[nnode - (nsteps-pivot_loc)/2 -1];
		}
		else
		{
			pn = &nodelist[nnode - (nsteps-pivot_loc-1)/2 -1];
			pn->Xe = p.dot(pn->Xe);
			pn->B.b = p.dot(pn->B.b);
		}

		while(pn != nodelist )
		{
			if(pn==nodelist) printf("DoPivot error\n");
			if (pn->father->lchild==pn) 
			{
				pn->father->p = p.dot(pn->father->p);
			}
			pn = pn->father;
			UpdateNode(pn);
		}
	}

	bool PivotAttempt()
	{
		Proposal prop(nsteps);
		DoPivot(prop.pivot_loc,prop.op);
		if(CheckPivot(prop))
		{
			DoPivot(prop.pivot_loc,prop.invop);
			return true;
		}
		else return false;
	}

	bool CheckPivot(Proposal prop)
	{
		node* pn;
		if (prop.pivot_loc % 2 == 1)
		{
			pn = &nodelist[nnode - (nsteps-prop.pivot_loc)/2 -1];
		}
		else
		{
		    if(prop.pivot_loc != 0) pn = &nodelist[nnode - (nsteps-prop.pivot_loc+1)/2 -1];
			else return false;
		}
		if (pn->father->rchild==pn)
		{
			pn = pn->father;
			while(pn->father->rchild==pn) {pn = pn->father;}
		}
		while(pn != nodelist)
		{
			if(pn==nodelist) printf("checkpivot error\n");
			if(pn->father->lchild==pn) 
			{
				if(Intersect(pn,pn+1,pn->Xe,pn->father->p)) return true;
				else pn = pn->father;
			}
			else 
			{
				if(Intersect(pn-1,pn,(pn-1)->Xe,pn->father->p)) return true;
				else pn = pn->father;
			}
		}
		return false;
	}

	bool Intersect(node* lc,node* rc,GPoint<double>& relXe,Matrix& p)
	{
		ncheck++;
		if(BoxIntersec(lc->B,rc->B,relXe,p)==false) return false;
		if(lc->isleaf && rc->isleaf) return true;
		if(lc->npoint >= rc->npoint)
		{
			/* Split the left SAW-tree; compare the SAW-trees which are closest together on the chain first, as
			they are the most likely to intersect. */
			if(Intersect(lc->rchild,rc,relXe - lc->lchild->Xe,p)) return true;
			else return(Intersect(lc->lchild,rc,relXe,p));
		}
		else
		{
			/* Split the right SAW-tree; compare the SAW-trees which are closest together on the chain first, as
			they are the most likely to intersect. */
			if(Intersect(lc,rc->lchild,relXe,p)) return true;
			else return(Intersect(lc,rc->rchild,relXe + rc->lchild->Xe,p.dot(rc->p)));
		}
	}

	void LineInit()
	{
		for (int i=nnode-1;i>=leafstartid;i--)
		{
			int indpoint = i-leafstartid;
			nodelist[i].Xe.point_assign(2,0,0);
			nodelist[i].B.a.assign(1,0,0,0.4,0);
			nodelist[i].B.b.assign(2,0,0,0.4,0);
		}
		for (int i=leafstartid-1;i>=0;i--)
		{
			nodelist[i].B = BoxMerge(nodelist[i].lchild->B,nodelist[i].rchild->B,nodelist[i].lchild->Xe,nodelist[i].p);
			nodelist[i].Xe = nodelist[i].rchild->Xe + nodelist[i].lchild->Xe;
		}
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

	void UpdateNode(node* pn)
	{
		pn->B = BoxMerge(pn->lchild->B,pn->rchild->B,pn->lchild->Xe,pn->p);
		pn->Xe = pn->p.dot(pn->rchild->Xe) + pn->lchild->Xe;
	}
	
	int ncheck;
	const char *data_fname, *final_name;
	node* nodelist;
	//node* ptemp;
	int nsteps;
	int deep;
	int nnode;
	int leafstartid;
};
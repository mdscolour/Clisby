#include "tree.h"


tree::tree(int tsteps):data_fname("data"),final_name("FinalWalk"),ncheck(0) //n is steps
{
	initialize(tsteps);
}

void tree::initialize( int tsteps )
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

	zeroid = 0;
	for (int i=0;i<deep;i++)
	{
		zeroid += twopowof(i);
	}

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

void tree::LineInit()
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

void tree::Record()
{
	FILE *fptr;
	// record the "data"
	double endnorm = nodelist[0].Xe.norm();
	fptr = fopen(data_fname, "a");
	fprintf(fptr,"%14.10f    %14.10f \n",123.0,endnorm*endnorm);
	fclose(fptr);
}

void tree::WriteDown()
{
	FILE *fptr;
	fptr = fopen(final_name, "w");

	int i;
	fprintf(fptr, "%d  %lf  \n", nsteps, -1.0);
	for (i = 0; i <= nsteps; i++) GetStepi(i).print(fptr);

	fclose(fptr);
}

Sphere tree::GetStepi( int i )
{
	Sphere t;
	node* pn;
	if (i % 2 == 1)
	{
		int pid = zeroid + (i-1)/2;
		if(pid>nnode-1) pid -= (nsteps+1)/2;
		pn = nodelist+pid;
		t = pn->B.b;
	}
	else
	{
		int pid = zeroid + (i)/2;
		if(pid>nnode-1) pid -= (nsteps+1)/2;
		pn = nodelist+pid;
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

int tree::run( int discard, int outer, int inner)
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

int tree::GoNSteps( int num )
{
	int accept = 0;
	for (int i=0;i<num;i++)
	{
		if(PivotAttempt()) accept++;
	}
	return accept;
}

bool tree::PivotAttempt()
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

inline void tree::DoPivot( int pivot_loc,Matrix p )
{
	node* pn;
	if (pivot_loc % 2 == 1)
	{
		int pid = zeroid + (pivot_loc-1)/2;
		if(pid>nnode-1) pid -= (nsteps+1)/2;
		pn = nodelist+pid;
	}
	else
	{
		int pid = zeroid + (pivot_loc)/2;
		if(pid>nnode-1) pid -= (nsteps+1)/2;
		pn = nodelist+pid;
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
		UpdateNodeData(pn);
	}
}

bool tree::CheckPivot( Proposal prop )
{
	node* pn;
	if (prop.pivot_loc % 2 == 1)
	{
		int pid = zeroid + (prop.pivot_loc-1)/2;
		if(pid>nnode-1) pid -= (nsteps+1)/2;
		pn = nodelist+pid;
	}
	else
	{
		if(prop.pivot_loc != 0) 
		{
			int pid = zeroid + (prop.pivot_loc)/2;
			if(pid>nnode-1) pid -= (nsteps+1)/2;
			pn = nodelist+pid-1;     // minus 1 here
		}
		else return false;
	}
	if (pn->father->rchild==pn)
	{
		pn = pn->father;
		while(pn->father->rchild==pn) {pn = pn->father;if(pn==nodelist) printf("checkpivot error 111\n");}
	}
	while(pn != nodelist)
	{
		if(pn==nodelist) printf("checkpivot error 222\n");
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

bool tree::Intersect( node* lc,node* rc,GPoint<double>& relXe,Matrix& p )
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

Matrix tree::Proposal::RandMatrix()
{
	return RefMatrix().dot(RMatrix());
}
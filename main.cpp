/***********************************************
so far this is nothing but a copy of the alg. in Clisby's paper
***********************************************/

#include "clis.h"

Point origin("zero");
Matrix I("identity");

void Merge(tree* wl,tree* wr,tree* w)
{
	w->n = wl->n+wr->n;
	w->B = wl->B.bound(wr->B.apply(w->q,wl->Xe));
	w->Xe = wl->Xe+w->q.dot(wr->Xe);
	w->X = wl->X+w->q.dot(wr->X)+wl->Xe*wr->n;
	w->X2 = wl->X2+wr->X2+wl->Xe.dot(w->q.dot(wr->X))*2+wl->Xe.dot(wl->Xe)*wr->n;
	return;
}

void LR(tree* w)
{
	tree* wt;
	Matrix qt;
	wt = w->wr;
	w->wr = wt->wr;
	wt->wr = wt->wl;
	wt->wl = w->wl;
	w->wl = wt;
	qt = w->q;
	w->q = qt.dot(w->wl->q);
	w->wl->q = qt;
	Merge(w->wl->wl,w->wl->wr,w->wl);
	return;
}
void RR(tree* w)
{
	tree* wt;
	Matrix qt;
	wt = w->wl;
	w->wl = wt->wl;
	wt->wl = wt->wr;
	wt->wr = w->wr;
	w->wr = wt;
	qt = w->q;
	w->q = w->wr->q;
	w->wr->q = w->q.inv().dot(qt);
	Merge(w->wr->wl,w->wr->wr,w->wr);
	return;
}

tree* Find_node(int nt,tree* w){while(w->wl!=0){w=w->wl;}return w+(nt-1);}
tree* Generate_SAW_tree(int n)
{
	int n_node = n,tolnum=n_node+n_node-1;
	int child_count=0,parent_count=n_node;
	tree* st = new tree[tolnum];
	for (int i=0;i<tolnum-1;i++)
	{
		st[i].q=I;
		if(i<n_node)
		{
			st[i].wl=0;
			st[i].wr=0;
			st[i].n=1;
			st[i].B=Box(1,1,0,0,0,0);
			st[i].Xe = Point(1,0,0);
			st[i].X = Point(1,0,0);
			st[i].X2 = 1;
		}
		else
		{
			st[i].wl=st+child_count;
			child_count++;
			st[i].wr=st+child_count;
			child_count++;
			Merge(st[i].wl,st[i].wr,st+i);
		}
		st[i].wp=st+parent_count;
		if(i%2==1)parent_count++;
	}
	st[tolnum-1].q=I;
	st[tolnum-1].wl=st+child_count;
	child_count++;
	st[tolnum-1].wr=st+child_count;
	child_count++;
	Merge(st[tolnum-1].wl,st[tolnum-1].wr,st+tolnum-1);
	st[tolnum-1].wp=0;
	return &st[tolnum-1];
}

int Random_integer_uniform(int a,int b){return a+int((b-a)*RNG_NAME());}
int Random_integer_log(int a,int b)//not efficiency
{
	if(b==a)return 0;
	int t,i;
	for (i=0;i<1000;i++)
	{
		t = a+int((b-a)*RNG_NAME());
		if(RNG_NAME()<=log(1+1.0/(t-a+1))/log(b-a+1.0)) break;
	}
	if(i==999)printf("error(Random_integer_log): 1000 trial not success but still in domain.\n");
	return t;
}

Matrix Random_symmetry()
{
	double ta,phi;
	ta = M_PI * RNG_NAME();
	phi = 2* M_PI * RNG_NAME();
	double l =sin(ta)*cos(phi),m = sin(ta)*sin(phi), n = cos(ta);

	double theta = 2* M_PI * RNG_NAME();
	Matrix rot(l*l*(1-cos(theta))+cos(theta), m*l*(1-cos(theta))-n*sin(theta), n*l*(1-cos(theta))+m*sin(theta), 
		l*m*(1-cos(theta))+n*sin(theta), m*m*(1-cos(theta))+cos(theta), n*m*(1-cos(theta))-l*sin(theta),
		l*n*(1-cos(theta))-m*sin(theta), m*n*(1-cos(theta))+l*sin(theta), n*n*(1-cos(theta))+cos(theta));

	ta = M_PI * RNG_NAME();
	phi = 2* M_PI * RNG_NAME();
	double a =sin(ta)*cos(phi),b = sin(ta)*sin(phi), c = cos(ta);
	Matrix ref(1-2*a*a, -2*a*b, -2*a*c, -2*a*b, 1-2*b*b, -2*b*c, -2*a*c, -2*b*c, 1-2*c*c);
	return ref.dot(rot);
}

void Accumulate_statistics(tree* w,bool success)
{
	//if(!success) return;
	FILE *fptr;
	char buffer[50]; // <- danger, only storage for 256 characters.
	sprintf_s(buffer, "data.txt");
	fptr = fopen(buffer, "a");

	fprintf(fptr,"%14.10f \n",w->Xe.norm());

	fclose(fptr);
}

void Shuffle_up(int n0,tree* w)
{
	if (n0==w->wl->n) return;
	else if(n0<w->wl->n)
	{
		Shuffle_up(n0,w->wl);
		RR(w);
	}
	else if(n0>w->wl->n)
	{
		Shuffle_up(n0 - w->wl->n,w->wr);
		LR(w);
	}
	return;
}
void Shuffle_down(tree* w)
{
	int nt = int((w->n+1)/2);
	if(nt == w->wl->n) return;
	else if (nt<w->wl->n)
	{
		RR(w);
		Shuffle_down(w->wr);
	}
	else if (nt>w->wl->n)
	{
		LR(w);
		Shuffle_down(w->wl);
	}
	return;
}
bool Intersect(Point xla,Matrix qla,tree* wl,Point xra,Matrix qra,tree* wr )
{
	Box blt = wl->B.apply(wl->q,xla);
	Box brt = wr->B.apply(wr->q,xra);
	if(blt.intersec(brt).notempty()==false) return false;
	if(wl->n<=2 && wr->n<=2) return true;
	if(wl->n>=wr->n)
	{
		if(Intersect(xla+qla.dot(wl->wl->Xe),qla.dot(wl->q),wl->wr,xra,qra,wr)) return true;
		return Intersect(xla,qla,wl->wl,xra,qra,wr);
	}
	else
	{
		if(Intersect(xla,qla,wl,xra,qra,wr->wl)) return true;
		return Intersect(xla,qla,wl,xra+qra.dot(wr->wl->Xe),qra.dot(wr->q),wr->wr);
	}
}

bool Shuffle_intersect(tree* w,Matrix q0,int wlc, int ilc)
{
	if(wlc==1){if(Intersect(origin,I,w->wl,w->wl->Xe+w->q.dot(q0.dot(w->wr->wl->Xe)),w->q.dot(q0.dot(w->wr->q)),w->wr->wr))return true;}
	else if(wlc==0){if(Intersect(origin,I,w->wl->wl,w->wl->Xe,w->q.dot(q0),w->wr))return true;}
	else if(wlc==-1){if(Intersect(origin,I,w->wl,w->wl->Xe,w->q.dot(q0),w->wr))return true;}

	if(w->wp == 0) return false;
	int ilc_new;
	if(w->wp->wp==0){}
	else if(w->wp->wp->wl==w->wp) ilc_new=true;
	else ilc_new = false;

	tree wt = *(w->wp);
	if(ilc) RR(&wt);
	else LR(&wt);

	return Shuffle_intersect(&wt,q0,ilc,ilc_new);
}

bool Attempt_pivot_simple(tree* w,int nt,Matrix qt)
{
	Shuffle_up(nt,w);
	w->q = w->q.dot(qt);
	bool intersection=Intersect(origin,I,w->wl,w->wl->Xe,w->q,w->wr);
	if(intersection) w->q = w->q.dot(qt.inv());
	else Merge(w->wl,w->wr,w);
	Shuffle_down(w);
	return (!intersection);
}

bool Attempt_pivot_fast(tree* w,int nt,Matrix qt)
{
	tree* wt = Find_node(nt,w);
	bool ilc;
	if(wt->wp->wl==0) {ilc=true;}
	else if(wt->wp->wl==wt){ilc=true;}
	else ilc=false;

	bool intersection = Shuffle_intersect(wt,qt,-1,ilc);
	if(!intersection)
	{
		Shuffle_up(nt,w);
		w->q = w->q.dot(qt);
		Shuffle_down(w);
		Merge(w->wl,w->wr,w);
	}
	return (!intersection);
}

void Pseudo_dimerize(tree* w)
{
	if(w->wl->full())Pseudo_dimerize(w->wl);
	if(w->wr->full())Pseudo_dimerize(w->wr);

	int nt;
	Matrix q,qt;
	int k;
	for(k=0;k<1000;k++)
	{
		if (w->wl->full())
		{
			nt = w->wl->n - Random_integer_log(1,w->wl->n);
			qt = Random_symmetry();
			Attempt_pivot_simple(w->wl,nt,qt);
		}
		if (w->wr->full())
		{
			nt = Random_integer_log(1,w->wr->n);
			qt = Random_symmetry();
			Attempt_pivot_simple(w->wr,nt,qt);
		}
		q = Random_symmetry();
		if(!Intersect(origin,I,w->wl,w->wl->Xe,w->q,w->wr)) break;
	}
	if(k==999)printf("error(Pseudo_dimerize): 1000 trial not success\n");

	Merge(w->wl,w->wr,w);

	for (int i=1;i<=sqrt(double(w->n));i++)
	{
		nt = w->wl->n - Random_integer_log(1,w->wl->n);
		qt = Random_symmetry();
		Attempt_pivot_simple(w,nt,qt);
		nt = w->wl->n + Random_integer_log(0,w->wr->n);
		qt = Random_symmetry();
		Attempt_pivot_simple(w,nt,qt);
	}
	return;
}

int main(int argc,char *argv[])
{
	int n = 1000;
	int n_dis=10,n_sample=10;
	tree *w = Generate_SAW_tree(n);
	Pseudo_dimerize(w);
	for(int i=1;i<=n_dis;i++)
	{
		int nt = Random_integer_uniform(1,n);
		Matrix qt = Random_symmetry();
		Attempt_pivot_simple(w,nt,qt);
		printf("%d\n",i);
	}
	for (int i=1;i<=n_sample;i++)
	{
		int nt = Random_integer_uniform(1,n);
		Matrix qt = Random_symmetry();
		bool success = Attempt_pivot_simple(w,nt,qt);
		Accumulate_statistics(w,success);
		printf("%d %d %d \n",i,nt,success);
	}

#ifdef _WIN32
	printf("\n");
	system("pause");
#endif
	return 0;
}
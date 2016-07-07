/***********************************************
so far this is nothing but a copy of the alg. in Clisby's paper
***********************************************/

#include "global.h"
#include "tree.h"

int main(int argc,char *argv[])
{
	tree A(1000-1);
	A.run(20000,1000,500);

#ifdef _WIN32
	printf("\n");
	system("pause");
#endif
	return 0;
}
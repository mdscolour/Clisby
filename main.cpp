/***********************************************
so far this is nothing but a copy of the alg. in Clisby's paper
***********************************************/

#include "global.h"
#include "tree.h"

int main(int argc,char *argv[])
{
	tree A(100-1);
	A.run();


#ifdef _WIN32
	printf("\n");
	system("pause");
#endif
	return 0;
}
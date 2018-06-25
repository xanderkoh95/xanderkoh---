#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
int main()
{

int i,j,N=9;
for(i=0;i<N-1;++i)
{
     for(j=0;j<N-1;++j)	
     {
     	if(i<(N+3)/2 && j<(N+3)/2)
     	{
		 printf("#\n");
	    }
	    else
	    {
	    	printf("(%d,%d)\n",i,j);
		}
	 }
}
return 0;
}

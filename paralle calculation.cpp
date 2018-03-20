#include<stdio.h>
#include<stdlib.h>
#include<omp.h>
#include<math.h>
#include<time.h>

int main()
{
	int i,j=0;
	clock_t t1, t2;
	
for(i=0;i<=10;++i)
	{
		j += i;
		printf("i=%d j=%d thread=%d\n",i,j,omp_get_thread_num());
	}
	printf("sum(1..10) = %d\n",j);
	j=0;
	#pragma omp parallel for reduction(+:j)//prevent wrong calculation when swapping the thread 
	for(i=0;i<=10;++i)
	{
		j += i;
		printf("i=%d j=%d thread=%d\n",i,j,omp_get_thread_num());
	}
	printf("sum(1..10) = %d\n",j);

	
	

return 0;	
	
}

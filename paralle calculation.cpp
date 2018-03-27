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
	//	printf("i=%d j=%d thread=%d\n",i,j,omp_get_thread_num());
	}
//	printf("sum(1..10) = %d\n",j);
	j=0;
	#pragma omp parallel for reduction(+:j)//prevent wrong calculation when swapping the thread 
	for(i=0;i<=10;++i)
	{
		#pragma omp atomic //prevent RAM to save the value twice ,generally slower than others
		j += i;
	//	printf("i=%d j=%d thread=%d\n",i,j,omp_get_thread_num());
	}
	//printf("sum(1..10) = %d\n",j);

int JG[10];	
//sum 0-9
j=0;
#pragma omp parallel num_threads(10)private(i)//paralle works prevent i is the independent variable not for share
{
	
	i = omp_get_thread_num();
	JG[i]=j;
   #pragma omp atomic //single scale variable  using reduction is more efficient 
	j+=i;
	}	
printf("thread number j = %d\n",j);



return 0;	
	
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<time.h>
int main()
{
	clock_t t1, t2,t_1,t_2,tt1,tt2;
	double a=1.234, b=2.456;
	int i, j, k, N=100000000;
	t1 = clock();
	for(i=0;i<N;++i)
	{
		a = a+b;
		b = a-b;
		a = a-b;
	}
	t2 = clock();
	tt1=t2-t1;
	printf("(+,-) loop time:%f\n",tt1/(double) CLOCKS_PER_SEC);
	printf("(a,b)=%f %f\n", a,b);
	
	t_1 = clock();
	for(i=0;i<N;++i)
	{
		a = a+b;
		b = a-b;
		a = a-b;
		a = a+b;
		b = a-b;
		a = a-b;
	}
	t_2 = clock();
	tt2=t_2-t_1;
	printf("(+,-) 6 times %f\n",(t_2-t_1)/(double) CLOCKS_PER_SEC);
	printf("(a,b)=%f %f\n", a,b);
	printf("Real time for plus minus %f\n",(tt2-tt1)/(double) CLOCKS_PER_SEC/3.0);
	t1 = clock();
	for(i=0;i<N;++i)
	{
		a = a*b;
		b = a/b;
		a = a/b;
	}
	t2 = clock();
	tt1=t2-t1;
	printf("(*,/) loop time:%f\n",tt1/(double) CLOCKS_PER_SEC);
	printf("(a,b)=%f %f\n", a,b);
	
	t_1 = clock();
	for(i=0;i<N;++i)
	{
		a = a*b;
		b = a/b;
		a = a/b;
		a = a*b;
		b = a/b;
		a = a/b;
	}
	t_2 = clock();
	tt2=t_2-t_1;
	printf("(*,/) 2 times %f\n",(t_2-t_1)/(double) CLOCKS_PER_SEC);
	printf("(a,b)=%f %f\n", a,b);
	printf("Real time for times and divided %f\n",(tt2-tt1)/(double) CLOCKS_PER_SEC/3.0);
	t1 = clock();
	for(i=0;i<N;++i)
	{
		a = sin(b);
		b = sin(a);
	}
	t2 = clock();
	tt1=t2-t1;
	printf("(sin) loop time:%f\n",tt1/(double) CLOCKS_PER_SEC);
	printf("(a,b)=%f %f\n", a,b);
	
	t_1 = clock();
	for(i=0;i<N;++i)
	{
	    a = sin(b);
		b = sin(a);
		a = sin(b);
		b = sin(a);
	}
	t_2 = clock();
	tt2=t_2-t_1;
	printf("(sin) 2 times %f\n",(t_2-t_1)/(double) CLOCKS_PER_SEC);
	printf("(a,b)=%f %f\n", a,b);
	printf(" Real time for calculating Sin %f\n",(tt2-tt1)/(double) CLOCKS_PER_SEC/3.0);

	return;
} 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
int Gauss(double **U,double **F,int N);
double Res(double **R,double **U,double **F,int N);
int main()
{
int i, j, k, p, N, Max_Steps = 100, M;
	double h, **U, **F, **R, r, r0, r1;
	time_t t1, t2;
	
	p = 3;
	N = 1 << p;
	M = (N-1)*(N-1);
	printf("N = %d\n",N);
	printf("Total Number of Unknowns = %d\n", M);
	h = M_PI/N;
	U = (double **) malloc( (N-1)*sizeof(double*) );
	U[0] = (double *) malloc( (N-1)*(N-1)*sizeof(double) );
	for(i=1;i<N-1;++i) U[i] = U[i-1] + (N-1);
	F = (double **) malloc( (N-1)*sizeof(double*) );
	F[0] = (double *) malloc( (N-1)*(N-1)*sizeof(double) );
	for(i=1;i<N-1;++i) F[i] = F[i-1] + (N-1);
	R = (double **) malloc( (N-1)*sizeof(double*) );
	R[0] = (double *) malloc( (N-1)*(N-1)*sizeof(double) );
	for(i=1;i<N-1;++i) R[i] = R[i-1] + (N-1);
	
	for(i=0;i<N-1;++i)
	{
		for(j=0;j<N-1;++j)
		{
			U[i][j] = 0.0;
			F[i][j] = 1.0*h*h;
		}
	}
	r = r0 = Res(R, U, F, N);
	k = 0;
	t1 = clock();
	while(r > r0*1e-6 && k < Max_Steps)
	{
		Gauss(U, F, N);
		r1 = Res(R, U, F, N);
		k ++;
		printf("%d:%e %f\n",k,r1,r1/r);
		r = r1;
	}
	t2 = clock();
	
	for(i=0;i<N-1;++i)
	{
	    for(j=0;j<N-1;++j)
	    {
	    	printf("(%d,%d)%f\n",i+1,j+1,U[i][j]);
		}
	}
	
	printf("Total time: %f\n", 1.0*(t2-t1)/CLOCKS_PER_SEC);
	return 0;	
	}
int Gauss(double **U,double **F,int N)
{
	int i, j;
	double s;
	
	for(i=0;i<N-1;++i) 
	{
		for(j=0;j<N-1;++j)
		{
			s = F[i][j];
			if(j-1>=0) s -= U[i][j-1];
			if(j+1< N-1) s -= U[i][j+1];
			if(i-1>=0) s -= U[i-1][j];
			if(i+1< N-1) s -= U[i+1][j];
			U[i][j] = s/(-4.0);
		}
	}
	for(i=0;i<N-1;++i) 
	{
		for(j=N-2;j>=0;j--)
		{
			s = F[i][j];
			if(j-1>=0) s -= U[i][j-1];
			if(j+1< N-1) s -= U[i][j+1];
			if(i-1>=0) s -= U[i-1][j];
			if(i+1< N-1) s -= U[i+1][j];
			U[i][j] = s/(-4.0);
		}
	}
	for(i=N-2;i>=0;i--) 
	{
		for(j=0;j<N-1;++j)
		{
			s = F[i][j];
			if(j-1>=0) s -= U[i][j-1];
			if(j+1< N-1) s -= U[i][j+1];
			if(i-1>=0) s -= U[i-1][j];
			if(i+1< N-1) s -= U[i+1][j];
			U[i][j] = s/(-4.0);
		}
	}
	for(i=N-2;i>=0;i--) 
	{
		for(j=N-2;j>=0;j--)
		{
			s = F[i][j];
			if(j-1>=0) s -= U[i][j-1];
			if(j+1< N-1) s -= U[i][j+1];
			if(i-1>=0) s -= U[i-1][j];
			if(i+1< N-1) s -= U[i+1][j];
			U[i][j] = s/(-4.0);
		}
	}
		
	
	for(i=0;i<(N-1)/2;++i)
	{
	    for(j=0;j<(N-1)/2;++j)
	    {
	    U[i][j]=0;
		}
	}
	
	return 0;
}
double Res(double **R, double **U, double **F, int N)
{
	int i, j;
	double s, r=0.0;
	
	for(i=0;i<N-1;++i) 
	{
		for(j=0;j<N-1;++j)
		{
     	if(i>(N-1)/2 && j>(N-1)/2)	
		{
			s = F[i][j]-U[i][j]*(-4);
			if(j-1>=0) s -= U[i][j-1];
			if(j+1< N-1) s -= U[i][j+1];
			if(i-1>=0) s -= U[i-1][j];
			if(i+1< N-1) s -= U[i+1][j];
			if(fabs(s) > r) r = fabs(s);
			R[i][j] = s;
		}	
		
		}
	}

	
	return r;
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
int main()
{
	clock_t t1, t2;				// variables for computing clocks 
	double **A, *x, *b, T1;// x=1234, a location in memory, x[0] the value at memory 1234
	int i, j, N=10000;
    
	
	for(N=10;N<=10;N*=2)
	{
		A = (double **) malloc( N * sizeof(double*) );//A[0],A[1];A[2] unlocated  //"double" memory 
		                                              //take N*N in RAM
		                                              // start from A[0]
		                                              //A[0][0],A[0][1], undefine yet 
		A[0] = (double *) malloc( N*N*sizeof(double));
		
		for(i=1;i<N;++i) A[i] = A[0]+ i*N;
		
		x = (double *) malloc( N * sizeof(double) );
		b = (double *) malloc( N * sizeof(double) );
		
		
		//A[1]=A[0]+N ---->A[2]=A[1]+N cannot to parallize
		#pragma omp parallel  // set a initial value on every thread
		srand(time(NULL));       //after that get the random number
		{
			#pragma omp parallel for	
		for(i=0;i<N;++i)//A[1]=A[0]+N,A[2]=A[0]+2N can be parallize
		{
			#pragma omp parallel for 
			for(j=0;j<N;++j)
			{
				A[i][j] = rand();
			}
			x[i] = rand();
		}
		t1 = clock();
		for(i=0;i<N;++i) 
		{
			b[i] = 0.0;
			for(j=0;j<N;++j)
			{
				b[i] += A[i][j]*x[j];
			}
		}
	}
		t2 = clock();
		T1 = (t2-t1)/(double) CLOCKS_PER_SEC;
		for(i=0;i<N;++i)
		{
			printf("%d\n",b[i]);
			//for(j=0;j<N;j++)
			{
	        //printf("%d",A[i][j]);
		    }
		  //printf("\n");
		}
		
		free(b);
		free(x);
		free(A[0]);
		free(A);	
	} 

	return 0;
} 

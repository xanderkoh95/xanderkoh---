#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
int main()
{
	clock_t t1, t2;				// variables for computing clocks 
	double **A, *x, *b,*c;// x=1234, a location in memory, x[0] the value at memory 1234
	int i,j,N;
	double T1,T2 ;
	
	for(N=5000;N<=5000;N*=2)
	{
		A = (double **) malloc( N * sizeof(double*) );//A[0],A[1];A[2] unlocated  //"double" memory 
		                                              //take N*N in RAM
		                                              // start from A[0]
		                                              //A[0][0],A[0][1], undefine yet 
		A[0] = (double *) malloc( N*N*sizeof(double));// given initial value
		#pragma omp parallel for
		for(i=1;i<N;++i) A[i] = A[0] + i*N;//A[1]=A[0]+N,A[2]=A[0]+2N can be parallize
		                                //A[1]=A[0]+N ---->A[2]=A[1]+N cannot be parallize
		
		x = (double *) malloc( N * sizeof(double) );
		b = (double *) malloc( N * sizeof(double) );
		c = (double *) malloc( N * sizeof(double) );
		
		
		#pragma omp parallel    
		{
		srand(time(NULL));       //after that get the random number
		#pragma omp parallel for private(i,j)
		for(i=0;i<N;++i)
		    {
			
			for(j=0;j<N;++j)
			  {
				A[i][j] = rand();
		
			   }
	      	x[i] = rand();
			}
		
	    }
	   // for(i=0;i<N;++i)// check the answer of matrix A
		 //{
		 //	for(j=0;j<N;++j){
			 
		 //	printf("%f\t",A[i][j]);
		 //	}
		//	 printf("\n");
		 //}
		t1 = clock(); 
		for(i=0;i<N;++i) 
		{
			b[i] = 0.0;
			for(j=0;j<N;++j)
			{
				b[i] += A[i][j]*x[j];
			}
		}
		t2 = clock();
		//for(i=0;i<N;++i) 
		//{
		//	printf("%f\n",b[i]); 
		//}
		T1 = (t2-t1)/(double) CLOCKS_PER_SEC;
	    printf(" matrix time vector before parallel is :%f\n",T1);
	    
		t1 = clock();
		#pragma omp parallel for private(i,j)
		for(i=0;i<N;++i) 
		{
			c[i] = 0.0;
			
			for(j=0;j<N;++j)
			{
				c[i] += A[i][j]*x[j];
			}
		}
		t2 = clock();
		T2 = (t2-t1)/(double) CLOCKS_PER_SEC;
	    printf(" matrix time vector after parallel is :%f\n",T2);
		printf(" time reduce  :%f\n",T1-T2);
		for (i=0;i<N;++i)
		printf("%f\n",fabs(b[i]-c[i]));
		printf("time reduce = %f  persent ",T2*100/T1);
		free(b);
		free(x);
		free(A[0]);
		free(A);	
	} 

	return 0;
} 

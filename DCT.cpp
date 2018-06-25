#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
int DCT(double *y_re, double *y_im, double *x_re, double *x_im, int N);
int main()
{
	int i,N;
	//printf("N=");
	//scanf("%d",&N);
	N=8;
	double *y_re, *y_im, *x_re, *x_im,T;
	y_re = (double *) malloc( N * sizeof(double));
	y_im = (double *) malloc( N * sizeof(double));
	x_re = (double *) malloc( N * sizeof(double));
    x_im = (double *) malloc( N * sizeof(double));
	
	for(i=0;i<N;++i)
	{
		x_re[i] = i+1;
		x_im[i] = 0.0;
	}

	DCT(y_re, y_im, x_re, x_im, N);
	for(i=0;i<N;++i)
	{
		printf("%f\n ", y_re[i]);
	}
	free(y_re);
	free(y_im);
	free(x_re);
	free(x_im);
	 
}
int DCT(double *y_re, double *y_im, double *x_re, double *x_im, int N)
{    if(N==2) 
	{
		// y, y[0] = x[0]+x[1], y[1] = x[0] - x[1]
		y_re[0] = x_re[0] + x_re[1];
		y_im[0] = x_im[0] + x_im[1];
		y_re[1] = x_re[0] - x_re[1]; 
		y_im[1] = x_im[0] - x_im[1];
	} else 
	{
		int k;
		double *y_even_re, *y_even_im, *y_odd_re, *y_odd_im;
		double *x_even_re, *x_even_im, *x_odd_re, *x_odd_im;
		double w_re, w_im, w_N_re, w_N_im, a, b, temp;
		y_even_re = (double *) malloc( N * sizeof(double));
		y_even_im = (double *) malloc( N * sizeof(double));
		x_even_re = (double *) malloc( N * sizeof(double));
		x_even_im = (double *) malloc( N * sizeof(double));
		y_odd_re = (double *) malloc( N * sizeof(double));
		y_odd_im = (double *) malloc( N * sizeof(double));
		x_odd_re = (double *) malloc( N * sizeof(double));
		x_odd_im = (double *) malloc( N * sizeof(double));
		for(k=0;k<N/2;++k)
		{
			x_even_re[k] = x_re[2*k];
			x_even_im[k] = x_im[2*k];
			x_odd_re[k]  = x_re[2*k+1];
			x_odd_im[k]  = x_im[2*k+1];
		}
		DCT(y_even_re, y_even_im, x_even_re, x_even_im, N/2);
		DCT(y_odd_re, y_odd_im, x_odd_re, x_odd_im, N/2);
	
		w_N_re =  cos(M_PI/(2*N));
		w_N_im = -sin(M_PI/(2*N));
		w_re   = 1.0;
		w_im   = 0.0; 
		for(k=0;k<N/2;++k)
		{
			a = w_re*y_odd_re[k] - w_im*y_odd_im[k];
			b = w_re*y_odd_im[k] + w_im*y_odd_re[k];
			y_re[k]     = y_even_re[k] + a;
			y_im[k]     = y_even_im[k] + b;
			y_re[N/2+k] = y_even_re[k] - a;
			y_im[N/2+k] = y_even_im[k] - b;
			temp = w_re;
			w_re = w_re*w_N_re - w_im*w_N_im;
			w_im = temp*w_N_im + w_im*w_N_re;
		}
		free(y_even_re);
		free(x_even_re);
		free(y_even_im);
		free(x_even_im);
		free(y_odd_re);
		free(y_odd_im);
		free(x_odd_re);
		free(x_odd_im);
	}	
 }

       
	


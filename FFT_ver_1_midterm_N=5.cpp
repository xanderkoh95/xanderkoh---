#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <time.h>

int Fast_Fourier_Transform(double *y_re, double *y_im, double *x_re, double *x_im, int N);
int main()
{
	int i;
	int N=5;
	double y_re[N], y_im[N], x_re[N], x_im[N];
	for(i=0;i<N;++i)
	{
		x_re[i] = i+1;
		x_im[i] = 0.0;
	}
	Fast_Fourier_Transform(y_re, y_im, x_re, x_im, N);
	for(i=0;i<N;++i)
	{
		printf("%f + %f i\n", y_re[i], y_im[i]);
	}
	
}
int Fast_Fourier_Transform(double *y_re, double *y_im, double *x_re, double *x_im, int N)
{
	double w1_re,w1_im,w2_re,w2_im,w3_re,w3_im,w4_re,w4_im,w6_re,w6_im,w8_re,w8_im,w9_re,w9_im,w12_re,w12_im,w16_re,w16_im;
	w1_re=0.3090169944;
	w1_im=0.9510565163;
	w2_re=-0.8090169944;
	w2_im=0.5877852523;
	w3_re=w2_re;
	w3_im=-w2_im;
	w4_re=w1_re;
	w4_im=-w1_im;
	w6_re=w1_re;
	w6_im=w1_im;
	w8_re=w2_re;
	w8_im=-w2_im;
	w9_re=w1_re;
	w9_im=-w1_im;
	w12_re=w2_re;
	w12_im=w2_im;
	w16_re=w1_re;
	w16_im=w1_im;
	if(N==5) 
	{	
		y_re[0] = x_re[0] + x_re[1] + x_re[2] + x_re[3] + x_re[4];
		y_im[0] = x_im[0] - x_im[1] - x_im[2] - x_im[3] - x_im[4];
		y_re[1] = x_re[0] + w1_re*x_re[1]-w1_im*x_im[1] + w2_re*x_re[2]-w2_im*x_im[2] + w3_re*x_re[3]-w3_im*x_im[3] + w4_re*x_re[4] - w4_im*x_im[4];
		y_im[1] = x_im[0] - w1_im*x_re[1]-w1_re*x_im[1] - w2_im*x_re[2]-w2_re*x_im[2] - w3_im*x_re[3]-w3_re*x_im[3] - w4_im*x_re[4] - w4_re*x_im[4];
		y_re[2] = x_re[0] + w2_re*x_re[1]-w2_im*x_im[1] + w4_re*x_re[2]-w4_im*x_im[2] + w6_re*x_re[3]-w6_im*x_im[3] + w8_re*x_re[4] - w8_im*x_im[4];
		y_im[2] = x_im[0] - w2_im*x_re[1]+w2_re*x_im[1] - w4_im*x_re[2]+w4_re*x_im[2] - w6_im*x_re[3]+w6_re*x_im[3] - w8_im*x_re[4] - w8_re*x_im[4];
		y_re[3] = x_re[0] + w3_re*x_re[1]-w3_im*x_im[1] + w6_re*x_re[2]-w6_im*x_im[2] + w9_re*x_re[3]-w9_im*x_im[3] + w12_re*x_re[4] - w12_im*x_im[4]; 
		y_im[3] = x_im[0] - w3_im*x_re[1]+w3_re*x_im[1] - w6_im*x_re[2]+w6_re*x_im[2] - w9_im*x_re[3]+w9_re*x_im[3] - w12_im*x_re[4] - w12_re*x_im[4];
		y_re[4] = x_re[0] + w4_re*x_re[1]-w4_im*x_im[1] + w8_re*x_re[2]-w8_im*x_im[2] + w12_re*x_re[3]-w12_im*x_im[3] + w16_re*x_re[4] - w16_im*x_im[4];
		y_im[4] = x_im[0] - w4_im*x_re[1]+w4_re*x_im[1] - w8_im*x_re[2]+w8_re*x_im[2] - w12_im*x_re[3]+w12_re*x_im[3] - w16_im*x_re[4] - w16_re*x_im[4];
		
	} 
	else 
	{
		//N=5
		int k;
		double *y_50_re, *y_50_im, *y_51_re, *y_51_im,*y_52_re,*y_52_im,*y_53_re,*y_53_im,*y_54_re,*y_54_im;
		double *x_50_re, *x_50_im, *x_51_re, *x_51_im,*x_52_re,*x_52_im,*x_53_re,*x_53_im,*x_54_re,*x_54_im;
		double w_re, w_im, w_N_re, w_N_im, a, b , c , d , temp ;
		y_50_re = (double *) malloc( N/5 * sizeof(double));
		y_50_im = (double *) malloc( N/5 * sizeof(double));
		y_51_re = (double *) malloc( N/5 * sizeof(double));
		y_51_im = (double *) malloc( N/5 * sizeof(double));
		y_52_re = (double *) malloc( N/5 * sizeof(double));
		y_52_im = (double *) malloc( N/5 * sizeof(double));
		y_53_re = (double *) malloc( N/5 * sizeof(double));
		y_53_im = (double *) malloc( N/5 * sizeof(double));
		y_54_re = (double *) malloc( N/5 * sizeof(double));
		y_54_im = (double *) malloc( N/5 * sizeof(double));
		x_50_re = (double *) malloc( N/5 * sizeof(double));
		x_50_im = (double *) malloc( N/5 * sizeof(double));
		x_51_re = (double *) malloc( N/5 * sizeof(double));
		x_51_im = (double *) malloc( N/5 * sizeof(double));
		x_52_re = (double *) malloc( N/5 * sizeof(double));
		x_52_im = (double *) malloc( N/5 * sizeof(double));
		x_53_re = (double *) malloc( N/5 * sizeof(double));
		x_53_im = (double *) malloc( N/5 * sizeof(double));
		x_54_re = (double *) malloc( N/5 * sizeof(double));
		x_54_im = (double *) malloc( N/5 * sizeof(double));		
	
		for(k=0;k<N/5;++k)
		{   //N=3
			x_50_re[k]=x_re[5*k];
	     	x_50_im[k]=x_im[5*k];
			x_51_re[k]=x_re[5*k+1];
			x_51_im[k]=x_im[5*k+1];
			x_52_re[k]=x_re[5*k+2];
			x_52_im[k]=x_im[5*k+2];
			x_53_re[k]=x_re[5*k+3];
			x_53_im[k]=x_im[5*k+3];
			x_54_re[k]=x_re[5*k+4];
			x_54_im[k]=x_im[5*k+4];
			
		}
		//N=3
		Fast_Fourier_Transform(y_50_re, y_50_im, x_50_re, x_50_im, N/5);
		Fast_Fourier_Transform(y_51_re, y_51_im, x_51_re, x_51_im, N/5);
		Fast_Fourier_Transform(y_52_re, y_52_im, x_52_re, x_52_im, N/5);
		Fast_Fourier_Transform(y_53_re, y_53_im, x_53_re, x_53_im, N/5);
		Fast_Fourier_Transform(y_54_re, y_54_im, x_54_re, x_54_im, N/5);

		//N=2: y_k = even_k + w_N^k odd_k = even_k + (a + bi)
		//N=3 : y_k=y_30[k]+ w_N^k(y_31[k])+w_N^2k(y_32[k]) 
		//         =y_30[k]+(w_re+iw_im)*(y_31_re+iy_31_im)+(w_re+iw_im)^2*( y_32_re+iy_32_im)
		w_N_re =  cos(2.0*M_PI/N);
		w_N_im = -sin(2.0*M_PI/N);
		w_re   = 1.0;     // initial value
		w_im   = 0.0; 
		
		
		for(k=0;k<N/5;++k)
		{
		
		y_re[k]=y_50_re[k] + y_re[1] + x_re[2] + x_re[3] + x_re[4];
		y_im[k]=
		y_re[N/5+k]=
		y_im[N/5+k]=
		y_re[2*N/5+k]=
		y_im[2*N/5+k]=
		y_im[3*N/5+k]=
		y_im[3*N/5+k]=
		y_im[4*N/5+k]=
		y_im[4*N/5+k]=
		
		
		
			
		
			temp = w_re;
			w_re = w_re*w_N_re - w_im*w_N_im;
			w_im = temp*w_N_im + w_im*w_N_re;
		}
		//N=3
		free(y_50_re);
		free(x_50_re);
		free(y_50_im);
		free(x_50_im);
		free(y_51_re);
		free(y_51_im);
		free(x_51_re);
		free(x_51_im);
		free(y_52_re);
		free(y_52_im);
		free(x_52_re);
		free(x_52_im);
		free(y_53_re);
		free(y_53_im);
		free(x_53_re);
		free(x_53_im);
		free(y_54_re);
		free(y_54_im);
		free(x_54_re);
		free(x_54_im);		
	}

}

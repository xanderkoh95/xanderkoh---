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
		//N=3
		int k;
		double *y_30_re, *y_30_im, *y_31_re, *y_31_im,*y_32_re,*y_32_im;
		double *x_30_re, *x_30_im, *x_31_re, *x_31_im,*x_32_re,*x_32_im;
		double w_re, w_im, w_N_re, w_N_im, a, b , c , d , temp ;
		y_30_re = (double *) malloc( N/3 * sizeof(double));
		y_30_im = (double *) malloc( N/3 * sizeof(double));
		y_31_re = (double *) malloc( N/3 * sizeof(double));
		y_31_im = (double *) malloc( N/3 * sizeof(double));
		y_32_re = (double *) malloc( N/3 * sizeof(double));
		y_32_im = (double *) malloc( N/3 * sizeof(double));
		x_30_re = (double *) malloc( N/3 * sizeof(double));
		x_30_im = (double *) malloc( N/3 * sizeof(double));
		x_31_re = (double *) malloc( N/3 * sizeof(double));
		x_31_im = (double *) malloc( N/3 * sizeof(double));
		x_32_re = (double *) malloc( N/3 * sizeof(double));
		x_32_im = (double *) malloc( N/3 * sizeof(double));
	
		for(k=0;k<N/3;++k)
		{   //N=3
			x_30_re[k]=x_re[3*k];
	     	x_30_im[k]=x_im[3*k];
			x_31_re[k]=x_re[3*k+1];
			x_31_im[k]=x_im[3*k+1];
			x_32_re[k]=x_re[3*k+2];
			x_32_im[k]=x_im[3*k+2];
			
		}
		//N=3
		Fast_Fourier_Transform(y_30_re, y_30_im, x_30_re, x_30_im, N/3);
		Fast_Fourier_Transform(y_31_re, y_31_im, x_31_re, x_31_im, N/3);
		Fast_Fourier_Transform(y_32_re, y_32_im, x_32_re, x_32_im, N/3);
		//N=2: y_k = even_k + w_N^k odd_k = even_k + (a + bi)
		//N=3 : y_k=y_30[k]+ w_N^k(y_31[k])+w_N^2k(y_32[k]) 
		//         =y_30[k]+(w_re+iw_im)*(y_31_re+iy_31_im)+(w_re+iw_im)^2*( y_32_re+iy_32_im)
		w_N_re =  cos(2.0*M_PI/N);
		w_N_im = -sin(2.0*M_PI/N);
		w_re   = 1.0;     // initial value
		w_im   = 0.0; 
		
		
		for(k=0;k<N/3;++k)
		{
			a = w_re*y_31_re[k] - w_im*y_31_im[k];//real part
			b = w_re*y_31_im[k] + w_im*y_31_re[k];//img part
			c = pow(w_re,2)*y_32_re[k]-2*w_re*w_im*y_32_im[k]-pow(w_im,2)*y_32_re[k];  //real part
			d = pow(w_re,2)*y_32_im[k]+2*w_re*w_im*y_32_re[k]-pow(w_im,2)*y_32_im[k];  
			
		/*	y_re[k]     = y_even_re[k] + a;
			y_im[k]     = y_even_im[k] + b;
			y_re[N/2+k] = y_even_re[k] - a;
			y_im[N/2+k] = y_even_im[k] - b;
			
		*/
		y_re[k]=y_30_re[k]+a+c;
		y_im[k]=y_30_im[k]+b+d;
		y_re[N/3+k]=y_30_re[k]-0.5*a+(pow(3.0,0.5)/2)*b-0.5*c-(pow(3.0,0.5)/2)*d;// problem in this one 
		y_im[N/3+k]=y_30_im[k]-0.5*b-(pow(3.0,0.5)/2)*a-0.5*d+(pow(3.0,0.5)/2)*c;// problem in this one
		y_re[2*N/3+k]=y_30_re[k]-0.5*a-(pow(3.0,0.5)/2)*b-0.5*c+(pow(3.0,0.5)/2)*d;// problem in this one
		y_im[2*N/3+k]=y_30_im[k]-0.5*b+(pow(3.0,0.5)/2)*a-0.5*d-(pow(3.0,0.5)/2)*c;// problem in this one
			
		
			temp = w_re;
			w_re = w_re*w_N_re - w_im*w_N_im;
			w_im = temp*w_N_im + w_im*w_N_re;
		}
		//N=3
		free(y_30_re);
		free(x_30_re);
		free(y_30_im);
		free(x_30_im);
		free(y_31_re);
		free(y_31_im);
		free(x_31_re);
		free(x_31_im);
		free(y_32_re);
		free(y_32_im);
		free(x_32_re);
		free(x_32_im);
		
	}

}

#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <time.h> 
int Fast_Fourier_Transform(double *y_re, double *y_im, double *x_re, double *x_im, int N);

int main()
{
	int i,N;
	N=43046721;
	double *y_re, *y_im, *x_re, *x_im,T;
	y_re = (double *) malloc( N * sizeof(double));
	y_im = (double *) malloc( N * sizeof(double));
	x_re = (double *) malloc( N * sizeof(double));
    x_im = (double *) malloc( N * sizeof(double));
	clock_t t1, t2;
	for(i=0;i<N;++i)
	{
		x_re[i] = i+1;
		x_im[i] = 0.0;
	}
	t1=clock();
	Fast_Fourier_Transform(y_re, y_im, x_re, x_im, N);
	t2=clock();
	T = 1.0*(t2-t1)/(double)CLOCKS_PER_SEC; 
	printf("%f s\n",T);
    for(i=0;i<N;++i)
	{
	//	printf("%f + %f i\n", y_re[i], y_im[i]);
	}
	free(y_re);
	free(y_im);
	free(x_re);
	free(x_im);
	 
}
int Fast_Fourier_Transform(double *y_re, double *y_im, double *x_re, double *x_im, int N)
{
	if(N%5 == 0)
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
		y_im[1] = x_im[0] - w1_im*x_re[1]+w1_re*x_im[1] - w2_im*x_re[2]+w2_re*x_im[2] - w3_im*x_re[3]+w3_re*x_im[3] - w4_im*x_re[4] + w4_re*x_im[4];
		y_re[2] = x_re[0] + w2_re*x_re[1]-w2_im*x_im[1] + w4_re*x_re[2]-w4_im*x_im[2] + w6_re*x_re[3]-w6_im*x_im[3] + w8_re*x_re[4] - w8_im*x_im[4];
		y_im[2] = x_im[0] - w2_im*x_re[1]+w2_re*x_im[1] - w4_im*x_re[2]+w4_re*x_im[2] - w6_im*x_re[3]+w6_re*x_im[3] - w8_im*x_re[4] + w8_re*x_im[4];
		y_re[3] = x_re[0] + w3_re*x_re[1]-w3_im*x_im[1] + w6_re*x_re[2]-w6_im*x_im[2] + w9_re*x_re[3]-w9_im*x_im[3] + w12_re*x_re[4] - w12_im*x_im[4]; 
		y_im[3] = x_im[0] - w3_im*x_re[1]+w3_re*x_im[1] - w6_im*x_re[2]+w6_re*x_im[2] - w9_im*x_re[3]+w9_re*x_im[3] - w12_im*x_re[4] + w12_re*x_im[4];
		y_re[4] = x_re[0] + w4_re*x_re[1]-w4_im*x_im[1] + w8_re*x_re[2]-w8_im*x_im[2] + w12_re*x_re[3]-w12_im*x_im[3] + w16_re*x_re[4] - w16_im*x_im[4];
		y_im[4] = x_im[0] - w4_im*x_re[1]+w4_re*x_im[1] - w8_im*x_re[2]+w8_re*x_im[2] - w12_im*x_re[3]+w12_re*x_im[3] - w16_im*x_re[4] + w16_re*x_im[4];
		
	} 
	else 
	{
		//N=5
		int k;
		double *y_50_re, *y_50_im, *y_51_re, *y_51_im,*y_52_re,*y_52_im,*y_53_re,*y_53_im,*y_54_re,*y_54_im;
		double *x_50_re, *x_50_im, *x_51_re, *x_51_im,*x_52_re,*x_52_im,*x_53_re,*x_53_im,*x_54_re,*x_54_im;
		double w_re, w_im, w_N_re, w_N_im, a, b , c , d , e , f , g , h , temp ;
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
		//N=5
		Fast_Fourier_Transform(y_50_re, y_50_im, x_50_re, x_50_im, N/5);
		Fast_Fourier_Transform(y_51_re, y_51_im, x_51_re, x_51_im, N/5);
		Fast_Fourier_Transform(y_52_re, y_52_im, x_52_re, x_52_im, N/5);
		Fast_Fourier_Transform(y_53_re, y_53_im, x_53_re, x_53_im, N/5);
		Fast_Fourier_Transform(y_54_re, y_54_im, x_54_re, x_54_im, N/5);

		//N=2: y_k = even_k + w_N^k odd_k = even_k + (a + bi)
		//N=3 : y_k=y_30[k]+ w_N^k(y_31[k])+w_N^2k(y_32[k]) 
		//         =y_30[k]+(w_re+iw_im)*(y_31_re+iy_31_im)+(w_re+iw_im)^2*( y_32_re+iy_32_im)
		//N=5 : y_k=y_50[k]+(w_re+iw_im)(y_51[k])+(w_re+iw_im)^2*( y_52_re+iy_52_im)+(w_re+iw_im)^3( y_53_re+iy_53_im)+(w_re+iw_im)^4( y_54_re+iy_54_im)
		w_N_re =  cos(2.0*M_PI/N);
		w_N_im = -sin(2.0*M_PI/N);
		w_re   = 1.0;     // initial value
		w_im   = 0.0; 
		
		
		for(k=0;k<N/5;++k)
		{
		a = w_re*y_51_re[k] - w_im*y_51_im[k];//real part
		b = w_re*y_51_im[k] + w_im*y_51_re[k];//img part
	    c = pow(w_re,2)*y_52_re[k]-2*w_re*w_im*y_52_im[k] - pow(w_im,2)*y_52_re[k];  //real part
		d = pow(w_re,2)*y_52_im[k]+2*w_re*w_im*y_52_re[k] - pow(w_im,2)*y_52_im[k]; 
		e = pow(w_re,3)*y_53_re[k]-3*w_re*pow(w_im,2)*y_53_re[k]-3*pow(w_re,2)*w_im*y_53_im[k]+pow(w_im,3)*y_53_im[k];
		f = pow(w_re,3)*y_53_im[k]-3*w_re*pow(w_im,2)*y_53_im[k]+3*pow(w_re,2)*w_im*y_53_re[k]-pow(w_im,3)*y_53_re[k];
		g = pow(w_re,4)*y_54_re[k]-4*pow(w_re,3)*w_im*y_54_im[k]-6*pow(w_re,2)*pow(w_im,2)*y_54_re[k]+4*w_re*pow(w_im,3)*y_54_im[k]+pow(w_im,4)*y_54_re[k];
		h = pow(w_re,4)*y_54_im[k]+4*pow(w_re,3)*w_im*y_54_re[k]-6*pow(w_re,2)*pow(w_im,2)*y_54_im[k]-4*w_re*pow(w_im,3)*y_54_re[k]+pow(w_im,4)*y_54_im[k];
		
		y_re[k]=y_50_re[k] + a + c + e + g ;
		y_im[k]=y_50_im[k] + b + d + f + h ;
		y_re[N/5+k]= y_50_re[k] + w1_re*a - w1_im*b + w2_re*c - w2_im*d + w3_re*e - w3_im*f + w4_re*g - w4_im*h; 
		y_im[N/5+k]= y_50_im[k] - w1_im*a + w1_re*b - w2_im*c + w2_re*d - w3_im*e + w3_re*f - w4_im*g + w4_re*h; 
     	y_re[2*N/5+k]= y_50_re[k] + w2_re*a - w2_im*b + w4_re*c - w4_im*d + w6_re*e - w6_im*f + w8_re*g - w8_im*h ;
		y_im[2*N/5+k]= y_50_im[k] - w2_im*a + w2_re*b - w4_im*c + w4_re*d - w6_im*e + w6_re*f - w8_im*g + w8_re*h;;
		y_re[3*N/5+k]= y_50_re[k] + w3_re*a - w3_im*b + w6_re*c - w6_im*d + w9_re*e - w9_im*f + w12_re*g - w12_im*h;
		y_im[3*N/5+k]= y_50_im[k] - w3_im*a + w3_re*b - w6_im*c + w6_re*d - w9_im*e + w9_re*f - w12_im*g + w12_re*h;
		y_re[4*N/5+k]= y_50_re[k] + w4_re*a - w4_im*b + w8_re*c - w8_im*d + w12_re*e - w12_im*f + w16_re*g - w16_im*h;
		y_im[4*N/5+k]= y_50_im[k] - w4_im*a + w4_re*b - w8_im*c + w8_re*d - w12_im*e + w12_re*f - w16_im*g + w16_re*h;
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
	
	

     else if (N%2 == 0)
     {
     if(N==2) 
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
		y_even_re = (double *) malloc( N/2 * sizeof(double));
		y_even_im = (double *) malloc( N/2 * sizeof(double));
		x_even_re = (double *) malloc( N/2 * sizeof(double));
		x_even_im = (double *) malloc( N/2 * sizeof(double));
		y_odd_re = (double *) malloc( N/2 * sizeof(double));
		y_odd_im = (double *) malloc( N/2 * sizeof(double));
		x_odd_re = (double *) malloc( N/2 * sizeof(double));
		x_odd_im = (double *) malloc( N/2 * sizeof(double));
		for(k=0;k<N/2;++k)
		{
			x_even_re[k] = x_re[2*k];
			x_even_im[k] = x_im[2*k];
			x_odd_re[k]  = x_re[2*k+1];
			x_odd_im[k]  = x_im[2*k+1];
		}
		Fast_Fourier_Transform(y_even_re, y_even_im, x_even_re, x_even_im, N/2);
		Fast_Fourier_Transform(y_odd_re, y_odd_im, x_odd_re, x_odd_im, N/2);
		// y_k = even_k + w_N^k odd_k = even_k + (a + bi)
		w_N_re =  cos(2.0*M_PI/N);
		w_N_im = -sin(2.0*M_PI/N);
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
       else if(N%3 ==0)
	   if(N==3) 
	{	
		y_re[0] = x_re[0] + x_re[1] + x_re[2];
		y_im[0] = x_im[0] + x_im[1] + x_im[2];
		y_re[1] = x_re[0] -0.5* x_re[1]+(pow(3.0,0.5)/2)* x_im[1]-0.5*x_re[2]-(pow(3.0,0.5)/2)*x_im[2]; 
		y_im[1] = x_im[0] -0.5* x_im[1]-(pow(3.0,0.5)/2)* x_re[1]-0.5*x_im[2]+(pow(3.0,0.5)/2)*x_re[2];
		y_re[2] = x_re[0] -0.5* x_re[1]-(pow(3.0,0.5)/2)* x_im[1]-0.5*x_re[2]+(pow(3.0,0.5)/2)*x_im[2];
		y_im[2] = x_im[0] -0.5* x_im[1]+(pow(3.0,0.5)/2)* x_re[1]-0.5*x_im[2]-(pow(3.0,0.5)/2)*x_re[2];
		
	} else 
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
		y_re[k]=y_30_re[k] + a + c;
		y_im[k]=y_30_im[k] + b + d;
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


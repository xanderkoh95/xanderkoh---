function u=dst_2(v)
n=length(v);
t=zeros(1,2*n+2);
t(2:1:n+1)=v;
t(n+3:1:2*n+2)=-flip(v);
b=fft(t);
c=imag(b)/2;
u=c(2:n+1);


end
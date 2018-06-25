function u=dct_2(v)

n=length(v);
t=zeros(1,4*n);
t(2:2:2*n)=v
g=fft(t);    
gg=real(g);
u=gg(1:n);

end
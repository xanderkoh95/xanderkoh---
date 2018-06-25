N=4;
number_U=(N-1)^2;
A=zeros(number_U);
Max_I=100000;
h=pi/(2*N);
for i=1:number_U
    for j=1:number_U
    if i==j
        A(i,j)=-4;
    elseif j==i+1 && mod(j,N-1)~=1;
        A(i,j)=1;
    elseif j==i-1 && mod(j,N-1)~=0;
        A(i,j)=1;
    elseif i==j-(N-1)
        A(i,j)=1;
    elseif i==j+(N-1)
        A(i,j)=1;
    end
end
end

D=diag(diag(A));
L=tril(A)-D;
P=triu(A)-D;
F=ones(number_U,1)*h^2;
U = zeros(number_U,1);
r0=norm(F);
for i=1:Max_I
U=D\(F-(P+L)*U);
if norm(F-A*U)/r0 < 1e-6
    i
    break;
    
end
end
mesh(reshape(U,N-1,N-1))



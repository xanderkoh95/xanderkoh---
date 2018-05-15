function[u]=CLT_A()%the size of the data which we sample (which is j (or N)here) increases
                   %the number of samples we take is fixed at 1000,
                   %create histogtams of mean_x for each J value (what do
                   %we see)?
for (j=10000:20000:100000)
    count=1;
    for k=1:1:1000
        N=j;
        x=1.75+0.5*rand(N,1);
        mean_X(count)=mean(x);
        var_(count)=var(x);
        count=count+1;
    end
end
plot (mean_X,j)
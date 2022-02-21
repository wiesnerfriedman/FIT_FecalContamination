function result=sedc4(D,a,zS,fC) %% sedc with flow connectivity

n=size(D,1); 
D2=-3*D;
Y=exp(D2./a).*fC;
kb=kron(ones(n,1),zS');
k=kb*Y';
result=sum(kron(k,2));
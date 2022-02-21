function result=sedc2(D,a,zS)

n=size(D,1); 
D2=-3*D;
Y=exp(D2./a);
kb=kron(ones(n,1),zS');
k=kb*Y';
result=sum(kron(k,2));
result=result';
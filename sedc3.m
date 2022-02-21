function result=sedc3(D,a)

%D2=-3*D;
% result=0.05.*exp(-3*D/(a))+0.95.*exp(-3*(D.^2)/((a*.001).^2));

%% Original
D2=-3*D;
result=exp(D2/a);
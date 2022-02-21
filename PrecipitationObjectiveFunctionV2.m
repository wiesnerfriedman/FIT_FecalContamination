function [ObjectFun,coef1,coef2,P1,P2,pval1,pval2,mdl]=PrecipitationObjectiveFunctionV2(Precip,T,a,b,Y)

tp=size(Precip,2);
NewP=Precip;
Dit=(1:tp);
P1=sum(NewP.*exp(-3*Dit./a),2);
P2=sum(NewP.*exp(-3*Dit./b),2);
P1P2=zscore(P1).*zscore(P2);
mdl=fitlm([zscore(P1) P1P2 T],Y);
coef1=mdl.Coefficients.Estimate(2);
coef2=mdl.Coefficients.Estimate(3);
pval1=mdl.Coefficients.pValue(2);
pval2=mdl.Coefficients.pValue(3);
ObjectFun=coef2-coef1;
if ~(coef1>0 && pval1<0.05)
% if (coef1<=0 || pval1>=0.05)
    mdl=fitlm([zscore(P1) T],Y);
    coef1=mdl.Coefficients.Estimate(2);
    coef2=mdl.Coefficients.Estimate(3);
    pval1=mdl.Coefficients.pValue(2);
    pval2=mdl.Coefficients.pValue(3);
    ObjectFun=1/(exp(coef1)+0.1);
end



%% 0- min MSE where p1 + p1p2
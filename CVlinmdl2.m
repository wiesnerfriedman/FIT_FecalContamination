function CVstruct=CVlinmdl2(tbl,nkfold)
% this is a k fold cross validation on the stepwise selection or "Test" stage

CVstruct=struct();

datatbl=tbl;
respons=datatbl(:,end); % extract response data
data=datatbl(:,1:end-1); % extract selected vars data
rng(1, 'twister');
cvIDX = crossvalind('Kfold',size(respons,1),nkfold);% create training and testing indices for the data
for i=1:nkfold
    idxb=(cvIDX==i); % test set
    idxa=~idxb; % train set
    tbla=tbl(idxa,:); 
    Xb=data(idxb,:); 
    yb=table2array(respons(idxb,:));
    mdla{i}=stepwiselm(tbla,'constant','Upper','linear','Criterion','aic');
    ybhat=predict(mdla{i},Xb);
    MSEcv(i)=(1/numel(yb)).*(sum((yb-ybhat).^2));
end
avMSE=mean(MSEcv);
varMSE=var(MSEcv);
CVstruct.cvModels=mdla;
CVstruct.avMSE=avMSE;
CVstruct.varMSE=var(MSEcv);
CVstruct.MSEcv=MSEcv;

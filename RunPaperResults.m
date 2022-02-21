%% Load data
pathname=pwd; %the path where the file ARGdataPub was saved
q=sprintf('%s%s',pathname,'/ARGdataPub.mat'); 
load(q)

%% Run Stepwise regression 

mdlBBsed=stepwiselm(tblBBsed ,'constant','Upper','linear','Criterion','aic');
mdlBBwat=stepwiselm(tblBBwat ,'constant','Upper','linear','Criterion','aic');
mdlHBsed=stepwiselm(tblHBsed ,'constant','Upper','linear','Criterion','aic');
mdlHBwat=stepwiselm(tblHBwat ,'constant','Upper','linear','Criterion','aic');

%% Run cross-validation 

cvBBsed=CVlinmdl2(tblBBsed,6);
cvBBwat=CVlinmdl2(tblBBwat,6);
cvHBsed=CVlinmdl2(tblHBsed,6);
cvHBwat=CVlinmdl2(tblHBwat,6);


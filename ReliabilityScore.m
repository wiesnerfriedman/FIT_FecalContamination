function RelInfo=ReliabilityScore(multisS, multizS, zR, sR, mt, flowtype, useRoadNet, assumcase,lb,ub, usepenalty, cvid)

%% INPUTS 

% multisS   {1 x c}   array with spatial source locations of size 1 by c, where c 
%                     is the number of candidate databases to be tried out             
% multizS   {1 x c}   with source scale information 
% zR        [1 x n]   the sample data
% sR        [2 x n]   the spatial locations for the sample data
% mt        {1, 2, 3} 1- Euclidean, 2- ORF, 3- GORF
% flowtype  {0, 1}    0- upstream cumulative riverstream lengths, 1-Strahler Order
% useRoadNet {0, 1}   0- Use Euclidean distance proxy, 1, use road network 
% assumcase  {0,1,2}  0- Inform based on fit, 1- inform based on positive regression
%                     coefficient, inform based on negative regression
%                     coefficient 
% lb         [1 x 1, 2, or 3]  lower bounds for the inform 
% ub         [1 x 1, 2, or 3]  upper bounds for the inform 
% cvid       k        number of folds
%% Output 

% RelInfo   struct   stores all the information from this! 

for icand=1:size(multisS,2)
    sS=multisS{icand}; 
    zS=multizS{icand}; 
    [D,zS,zRq,~,~,modeltype]=KewauneeObjFunLoad(sS,zS,sR,zR,flowtype,mt,useRoadNet);
    rng(1,'twister'); % ensure that selection occurs pseudo randomly
    indices=crossvalind('Kfold',zR,cvid);
    for cvi=1:cvid
        qcvid=(indices~=cvi);
        notqcvid=(indices==cvi);
        params{cvi}=SEDCparamOptim(D,zS,zRq,0,assumcase,lb,ub,mt,qcvid,0,usepenalty); % get the hyperparameters
        [~,Coef(cvi)]=KewauneeObjectiveFunction1(D,zS,params{cvi}(1,1),params{cvi}(1,2),params{cvi}(1,3),0,zRq,assumcase,notqcvid,modeltype,usepenalty);
    end
    SS(icand)=(sum(Coef>0)+1).*sum(Coef);
    ParamSD(icand)=sqrt(var(Coef));
    RS(icand)=SS(icand)./ ParamSD(icand);
    Coefficients{icand}=Coef; 
    Hyperparameters{icand}=params;
end
cands=1:size(multisS,2); 
Candidates=cands'; 
SignStability=SS'; 
HyperparameterSD=ParamSD';
Reliability=RS'; 


RelTab=table(Candidates, SignStability, HyperparameterSD, Reliability);
RelInfo=struct(); 

RelInfo.RelTab=RelTab; 
RelInfo.KFoldCoeffs=Coefficients; 
RelInfo.KFoldParams=Hyperparameters; 


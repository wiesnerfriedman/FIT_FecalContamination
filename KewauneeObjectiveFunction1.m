function [ObjectFun,Coef,pValue,R2,CIr,result,info,penalty]=KewauneeObjectiveFunction1(D,zS,a,b,c,STVs,yqi,assumcase,idxrun,modeltype,usepenalty)

%% this objective function is used in evaluating the optimal alpha value (s) for model types 1-5 of SEDCTOR


% INPUT

% D|modeltype=1or2[nxm]      n by m matrix of distances between each
%                           observation point and each source point
% D| modeltype=3-5{q}[nxm] qxr array of up to three n by m distance matrices
%                           and one n by m flow connectivity matrix, the
%                           flow connectivity matrix always comes last
% zS|modeltype=1  [nx1]      vector of size m of values associated with each
%                           source point
% a|modeltype=1or2scalar     scalar value representing the exponential
%                           distance decay range
% a|modeltype=3-5 vector     vector of values representing the exponential
%                           distance decay ranges for different processes
%                           higherarchy is 1)T 2)O 3)R
% STVs            [mxp]      matrix of selected temporal variables
% y|modeltype=1   [mx1]      vector of size m of observation values
% y|modeltype=2-5 [mx2]      vector of size m of observation values the flow at that point, Qi
% assumcase       {0,1,2}    value indicating assumptions about coefficien
%                           0- no assumptions 1- positive 2- negative
% idxrun          [mx1]      logical vector of values indicating which if
%                           just a 2, then it means no indexing

% subset data based on index idxrun

if numel(idxrun)>0
    idxrun=logical(idxrun);
    if STVs~=0
       STVs=STVs(idxrun,:);
    end 
    
    
    switch modeltype
        case 1 % Euclidean Distance SEDC-E (Euclidean)
            D=D{1}(idxrun,:);
            y=yqi(idxrun,1);
            result=sedc2(D,a,zS);
            penalty=0; 
        case 2 % Can use, but not used in manuscripts because case 3 is better
            y=yqi(idxrun,1);
            DIJR=D{1}(idxrun,:);
            fC=D{2}(idxrun,:);
            Qi=yqi(idxrun,2);
            Qi(Qi==0)=min(Qi(Qi~=0))/10000000;
            xc=sedc4(DIJR,a,zS,fC);
            result=xc./Qi';
            result=result';
            penalty=0;
        case 3 % Overland Flow and River Distance Flow Connected (ORF)
            y=yqi(idxrun,1);
            DIJO=D{1};
            DIJR=D{2}(idxrun,:);
            fC=D{3}(idxrun,:);
            Qi=yqi(idxrun,2);
            Qi(Qi==0)=min(Qi(Qi~=0))/10000000;
            A=sedc3(DIJO,a);
            B=sedc3(DIJR,b);
            xa=A.*B;
            xb=zS'.*xa.*fC;
            xc=sum(xb,2);
           % result=xc./Qi;
		    result=xc;
            penalty=0;   
        case 4 % Do not use
            y=yqi(idxrun,1);
            DGJ4=D{1};
            DIGR4=D{2}(idxrun,:);
            fC=D{3}(idxrun,:);
            Qi=yqi(idxrun,2);
            Qi(Qi==0)=min(Qi(Qi~=0))/10000000;
            A1=sedc3(DGJ4,a);
            A2=sum(A1,1);
            A2(A2==0)=1;
            A=sum((A1./A2).*zS',2);
            B=sedc3(DIGR4,b);
            C=sum(A'.*B.*fC,2);
            result=C./Qi;
        case 5 % Ground transport to location, overland flow, and river distance flow connected (GORF)
            if iscell(zS)
                y=yqi(:,1);
                %DGJ5=D{1};
                DGJ5=D{1};
                DIGO5=D{2};
                DIGR5=D{3};
                fC=D{4};
                Qi=yqi(:,2);
                Qi(Qi==0)=min(Qi(Qi~=0))/10000000;
                A=GORFtinypartv2(DGJ5,zS,a,1); 
                B=sedc3(DIGO5,b);
                C=sedc3(DIGR5,c);
                xa=sum(A.*B.*C.*fC,2);
                result=xa./Qi;
%                 A1=sedc3(DGJ5,a);
%                 A2=sum(A1.*zS{2},1);
%                 A2(A2==0)=1;
%                 A=sum(((A1.*zS{2})./A2).*zS{1}',2);
%                 B=sedc3(DIGO5,b);
%                 C=sedc3(DIGR5,c);
%                 xa=sum((A.*B)'.*(C.*fC),2);
%                 result=xa./Qi;
            else
                y=yqi(:,1);
                DGJ5=D{1};
                DIGO5=D{2};
                DIGR5=D{3};
                fC=D{4};
                Qi=yqi(:,2);
                Qi(Qi==0)=min(Qi(Qi~=0))/10000000;
                A=GORFtinypartv2(DGJ5,zS,a,0); 
                B=sedc3(DIGO5,b);
                C=sedc3(DIGR5,c);
                xa=sum((A.*B)'.*(C.*fC),2);
                result=xa./Qi;
                %result=xa;
            end
    end
else
    switch modeltype
        case 1 % Euclidean 
            y=yqi(:,1);
            D=D{1};
            result=sedc2(D,a,zS);
            penalty=0;  
        case 2 % RF
            y=yqi(:,1);
            DIJR=D{1};
            fC=D{2};
            Qi=yqi(:,2);
            Qi(Qi==0)=min(Qi(Qi~=0))/10000000;
            xc=sedc4(DIJR,a,zS,fC);
            result=xc./Qi';
            result=result';
            penalty=0;
        case 3 % ORF
            y=yqi(:,1);
            DIJO=D{1};
            DIJR=D{2};
            fC=D{3};
            Qi=yqi(:,2);
            Qi(Qi==0)=min(Qi(Qi~=0))/10000000;
            A=sedc3(DIJO,a);
            B=sedc3(DIJR,b);
            xa=A.*B;
            xb=zS'.*xa.*fC;
            xc=sum(xb,2);
            result=xc./Qi;
            %result=xc;
            penalty=0; 
        case 4 % Do not use
            y=yqi(:,1);
            DGJ4=D{1};
            DIGR4=D{2};
            fC=D{3};
            Qi=yqi(:,2);
            Qi(Qi==0)=min(Qi(Qi~=0))/10000000;
            A1=sedc3(DGJ4,a);
            A2=sum(A1,1);
            A2(A2==0)=1;
            A=sum((A1./A2).*zS',2);
            B=sedc3(DIGR4,b);
            C=sum(A'.*B.*fC,2);
            result=C./Qi;
            penalty=0;
        case 5 % GORF
            if iscell(zS)
                y=yqi(:,1);
                %DGJ5=D{1};
                DGJ5=D{1};
                DIGO5=D{2};
                DIGR5=D{3};
                fC=D{4};
                Qi=yqi(:,2);
                Qi(Qi==0)=min(Qi(Qi~=0))/10000000;
                A=GORFtinypartv2(DGJ5,zS,a,1); 
                B=sedc3(DIGO5,b);
                C=sedc3(DIGR5,c);
                xa=sum(A'.*B'.*C.*fC,2);
                result=xa./Qi;
%                 A1=sedc3(DGJ5,a);
%                 A2=sum(A1.*zS{2},1);
%                 A2(A2==0)=1;
%                 A=sum(((A1.*zS{2})./A2).*zS{1}',2);
%                 B=sedc3(DIGO5,b);
%                 C=sedc3(DIGR5,c);
%                 xa=sum((A.*B)'.*(C.*fC),2);
%                 result=xa./Qi;
            else
                y=yqi(:,1);
                DGJ5=D{1};
                DIGO5=D{2};
                DIGR5=D{3};
                fC=D{4};
                Qi=yqi(:,2);
                Qi(Qi==0)=min(Qi(Qi~=0))/10000000;
                A=GORFtinypartv2(DGJ5,zS,a,0); 
                B=sedc3(DIGO5,b);
                C=sedc3(DIGR5,c);
                xa=sum((A.*B)'.*(C.*fC),2);
                result=xa./Qi;
                %result=xa; 
            end
    end
end 
% verify that the result is not comprised of only zeros or only ones (not
% rank deficient)
str1=ones(size(result));
idx1=strcmp(char(result),char(str1)); % indicates that the column is filled with ones
%% try square root 
tr=result;
idxinf=isinf(tr);
idxnan=isnan(tr); 
tr(idxinf)=[];
if numel(tr)==0
    tr(idxnan)=[];
end

if size(STVs,1)>1
    STVs2=STVs(~idxnan,:);
end
tr2=result;
tr2(idxinf)=nanmax(tr); % verify there aren't inf or NaN issues
tr2(idxnan)=nanmean(tr);

% zscore the design matrix
tr2(~idxinf|~idxnan)=zscore(tr2(~idxinf|~idxnan)); 
if STVs~=0
    STVs2=zscore(STVs2);
    
    
    stvs3=STVs;
    stvs3((idxinf|idxnan),:)=[];
    fulmat=[tr stvs3];
else
    fulmat=tr;
end
%matrankidx=(rank(fulmat)==size(fulmat,2)); % verify again that the matrix is not rank deficient due to other causes

if (sum(result)==0) || (idx1==1) || (numel(unique(result))==1)  || sum(idxinf)==size(result,1) || sum(idxnan)==size(result,1) %|| (matrankidx==0) 
    switch assumcase
        case 0
            mdl=[];
            Coef=0;
            pValue=1;
            R2=0;
            MSE=1;
            CIs=[];
            CIr=[-1 1];
        case 1 % coefficient must be positive
            mdl=[];
            Coef1=-500000;
	    Coef=Coef1;
            pValue=1;
            R2=0;
            MSE=1;
            CIs=[];
            CIr=[-1 1];
	    vr=NaN; 
	    penalty=NaN; 
        case 2 % coefficient must be negative
            mdl=[];
            Coef1=500000;
	    Coef=Coef1;
            pValue=1;
            R2=0;
            MSE=1;
            CIs=[];
            CIr=[-1 1];
    end
else
    if size(STVs,1)>1
        X=[tr2 STVs2]; % see above -- standardized the variables using zscore, while ignoring inf
    else
        X=tr2;
    end
  %  mdl=fitlm(X,y);
  %  Coef=mdl.Coefficients.Estimate(2);
  %  pValue=mdl.Coefficients.pValue(2);
  %  R2=mdl.Rsquared.Ordinary;
  %  MSE=mdl.MSE;
  %  CIs=coefCI(mdl);
  %  CIr=CIs(2,:);

    uvals=unique(tr2,'rows');
    
    for iu=1:size(uvals,1)
	idx=(tr2==uvals(iu,:)); 
	nsamps=sum(idx(:,1)); 
	weighting(idx(:,1))=1./nsamps;
    end

    weighting=weighting./(sum(weighting))./size(y,1); 
    
    mdl=fitglm(X,y,'weights',weighting); 
    Coef=mdl.Coefficients.Estimate(2);
    pValue=mdl.Coefficients.pValue(2);
    R2=mdl.Rsquared.Ordinary;
    MSE=mdl.LogLikelihood;
    CIs=coefCI(mdl);
    CIr=CIs(2,:);
    penalty=0; 
    if usepenalty==1
    	if modeltype==5
        	penalty=KewauneePenaltyFunction1(a,b,2,100,8);
    	elseif modeltype==3
        	penalty=KewauneePenaltyFunction2(a,2,5,4);
    	end  
    end
%rng(0,'Twister') 
% bootstrap approach! 
%     mdlboot=BootStrapRegression(X,y,500,size(y,1));
%     Coef=mdlboot.Coefficients(2); 
%     pValue=mdlboot.pValue(2); 
%     R2=mdlboot.R2; 
%CIs=0;% if successful go back and change these.  
%CIr=0;

end
if abs(Coef)>100
Coef1=Coef./100000;
else
Coef1=Coef; 
end  
switch assumcase
    case 0
        ObjectFun=pValue; % do nothing; % reward results with high R2/pValue
    case 1 % coefficient must be positive
        ObjectFun=(1/exp(Coef1))+penalty;
    case 2 % coefficient must be negative
        ObjectFun=exp(Coef1);%10.^(Coef);
end


info=struct();
info.n=size(result,1); 
info.NaNvalues=sum(idxnan); 
info.INFvalues=sum(idxinf); 
%info.MSE=MSE;
%info.mdl=mdl; 
info.params=[a b c]; 


function [result,result2,result3]=KewauneeSEDCresults(D,zS,params,yqi,modeltype)

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

switch modeltype
    case 1 % Euclidean Distance SEDC-E
        a=params(1);
        D=D{1};
        y=yqi(:,1);
        result=sedc2(D,a,zS);
        result2=0;
        result3=0;
    case 2 % River Distance-Flow Connected SEDC-R
        a=params(1);
        y=yqi(:,1);
        DIJR=D{1};
        fC=D{2};
        Qi=yqi(:,2);
        Qi(Qi==0)=min(Qi(Qi~=0))/10000000;
        xc=sedc4(DIJR,a,zS,fC);
        result=xc./Qi';
        result=result';
        result2=0;
        result3=0;
    case 3 % Overland Flow and River Distance Flow Connected SEDC-OR
        a=params(1);
        b=params(2);
        y=yqi(:,1);
        DIJO=D{1};
        DIJR=D{2};
        fC=D{3};
        Qi=yqi(:,2);
        Qi(Qi==0)=min(Qi(Qi~=0))/10000000;
        A=sedc3(DIJO,a);
        B=sedc3(DIJR,b);
        if size(A,2)~=size(B,2)
            A=A';
        end 
        xa=A.*B;
        xb=zS'.*xa.*fC;
        xc=sum(xb,2);
        result=xc./Qi;
        result2=zS'.*A';
        result3=0;
    case 4 % Ground Transport to River and River Distance Flow Connected SEDC-TR
        a=params(1);
        b=params(2);
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
        result2=A;
        result3=0;
    case 5 % Ground transport to location, overland flow, and river distance flow connected SEDC-TOR
        a=params(1);
        b=params(2);
        c=params(3);
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
            result2=A;
            result3=(A.*B)'; 
        else
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
            result2=A;
            result3=(A.*B)'; 
        end
end


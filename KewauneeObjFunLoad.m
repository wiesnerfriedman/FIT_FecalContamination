function [D,zS,yqi,ub,lb,modeltype]=KewauneeObjFunLoad(sS,zS,sR,zR,flowtype,mt,useRoadNet)

load('KewauneeRiverNetwork_Minifile.mat','FI','fC','rD','SOflow','URCRM')

%URCRM=URCRM./max(URCRM); 

SOflow=SOflow; 
if flowtype==2
    Qall=SOflow;
elseif flowtype==1
    Qall=URCRM; 
end

RNpoints=FI(:,1:2); 
fCorig=fC; 
idxsK=size(sR,1)~=size(RNpoints,1);
if idxsK
    IDXsR=FindSitesonNetwork(RNpoints,sR,Qall);
    %IDXsR=knnsearch(RNpoints,sR); 
    if flowtype==2
        Qi=SOflow(IDXsR,1);
        yqi=[zR Qi];
    elseif flowtype==1
        Qi=URCRM(IDXsR,1);
        yqi=[zR Qi];
    end
else
    if flowtype==2
        Qi=SOflow;
        yqi=[zR Qi];
    elseif flowtype==1
        Qi=URCRM;
        yqi=[zR Qi];
    end
end

% Euclidean
Ebounds=[50 34000];
% River Distance
Rbounds=[249 1300000];
% Overland Flow
Obounds=[1 10000];
% Ground Transport
Gbounds=[1 340000];

switch mt
    case 1 % Euclidean
        DIJE=sectionedcoord2dist(sR,sS); %% sources to sampling point
        D{1}=DIJE;
        clear DIJE
        lb = Ebounds(1);
        ub = Ebounds(2);
        modeltype=1;
    case 2 % ORF
        DGJ4=sectionedcoord2dist(RNpoints,sS); %% source points to the river network
        [DIJO, IdxSS2RN]=min(DGJ4); %% overland distance to the nearest river network point and index indicating the rn point
%         [QID,QD]=knnsearch(RNpoints,sS,'k',3); 
%         IdxSS2RN=QID(:);
%         DIJO=QD(:);
        if ~idxsK
            DIJR=rD(:,IdxSS2RN); 
            fC=fCorig(:,IdxSS2RN);
        else 
            DIJR=rD(IDXsR,IdxSS2RN); %% from overland flow point from source to the river
            fC=fCorig(IDXsR,IdxSS2RN);
        end 
        D{1}=DIJO;
        D{2}=DIJR;
        D{3}=fC;
        clear fC DIJR DIJO IdxSS2RN DGJ4
        lb = [Obounds(1),Rbounds(1)];
        ub = [Obounds(2),Rbounds(2)];
        modeltype=3;
    case 3 % GORF
        sS1=sS{1}; % location initial
        sS2=sS{2}; % location transported to
        if useRoadNet==1
            DGJ5=KewRoadDistanceMatrix(sS); % distance between the source and the grid point (transport distance
        else
            DGJ5=sectionedcoord2dist(sS2,sS1); % distance between the source and the grid point (transport distance
        end
        [IdxSS2RN2,DIGO5]=knnsearch(RNpoints,sS2); % distance between the grid point and the nearest neighbor river network point with index to RN point
%         [QID,QD]=knnsearch(RNpoints,sS2,'k',3); 
%         tmp1=sS2(:,1); 
%         tmp2=sS2(:,2); 
%         q1=repmat(tmp1,[1 3]); 
%         q2=repmat(tmp2,[1 3]);
% %         sS2b=[q1(:) q2(:)];
%         IdxSS2RN2=QID(:);
%         DIGO5=QD(:);
%       DGJ5=sectionedcoord2dist(sS2b,sS1); % distance between the source and the grid point (transport distance
        if ~idxsK
            DIGR5=rD(:,IdxSS2RN2);
            fC5=fCorig(:,IdxSS2RN2);
        else
            DIGR5=rD(IDXsR,IdxSS2RN2);
            fC5=fCorig(IDXsR,IdxSS2RN2);
        end
        % set lowerbound and upper bounds
        lb = [Gbounds(1),Obounds(1),Rbounds(1)];
        ub = [Gbounds(2),Obounds(2),Rbounds(2)];
        % set the distance matrix
        D{1}=DGJ5;
        D{2}=DIGO5;
        D{3}=DIGR5;
        D{4}=fC5;
        clear DGJ5 IdxSS2RN2 DIGO5 DIGR5 fC5
        modeltype=5;
end
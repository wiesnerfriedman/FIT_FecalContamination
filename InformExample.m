%% Load data
pathname=pwd; %the path where the file ARGdataPub was saved
q=sprintf('%s%s',pathname,'/ARGdataPub.mat'); 
load(q)

%% ORF example set source and response data
% the source is high intensity developed land cover (HID) represented by 
% by equally weighted gridded points that intersect with the land cover 
% and the response is the log10 relative abundance of human bacteroides 
% in sediment 

% The spatial locations of the source
sS1=ExampleSources.HIDoptionA(:,1:2); 

% The scale of the source (here each equal to 1 because of being equally
% weighted gridded points)
zS1=ExampleSources.HIDoptionA(:,3); 

% The spatial locations of the response (sampling sites) 
sR1=ExampleResponses.RAHuBacSed(:,1:2);

% The response value 
zR1=ExampleResponses.RAHuBacSed(:,3);

%% load information needed for the inform 
% The model type needs to be set to 2
mt=2; 

% The type of proxy for flow we are using is Strahler Order flow, rather
% than upstream contributing river length
flowtype=2;

% We are not using the road network, instead using Euclidean distance proxy
useRoadNet=0;

% load the distance matrices and other needed information for the inform 
[D1,zS1,zRq1,ub1,lb1,modeltype1]=KewauneeObjFunLoad(sS1,zS1,sR1,zR1,flowtype,mt,useRoadNet);

%% get the parameter values 
% we inform independently of the meteorological factors 
STVs=0;

% we will use positive coefficients rather than R2 to inform
assumcase=1;

% we will set lower and upper bounds for the optimization based on the
% plots that were manually generated in the SI of the manuscript see 
% Figure S8(g) where we see more global optimums

% River Distance hyperparameter upper and lower bounds
Rbounds=[249 3000]; 

% Overland Flow hyper parameter upper and lower bounds
Obounds=[100 3000];

% lower bounds
lb = [Obounds(1),Rbounds(1)];
ub = [Obounds(2),Rbounds(2)];

% We don't want to subset the data in any way
idxset=[]; 

% We want to use just pattern search to get the parameters, not bruteforce
bruteforce=0; 

% We want to use the penalty function to direct the parameter values away
% from things that have poor regression qualities 
usepenalty=1;

% INFORM!  
params1=SEDCparamOptim(D1,zS1,zRq1,STVs,assumcase,lb,ub,mt,idxset,bruteforce,usepenalty);

% We view the parameters that we found 

sprintf('alpha_O=%g km, alpha_R=%g km', params1(1)/1000, params1(2)/1000)

%% obtain the HID source term values pre-z-scored

ORFhid=KewauneeSEDCresults(D1,zS1,params1,zRq1,modeltype1);

%%  GORF example set source and response data
% the source is AFO via hauling manure to application fields. AFOs are
% represented by WPDES CAFOs weighted by animal units. The location, rather
% than size of the fields matter for the manure application variable

% The spatial locations of the source
sS2{1}=ExampleSources.AFOoptionA(:,1:2); 

% The scale of the source (here as the animal units)
zS2=ExampleSources.AFOoptionA(:,3); 

% The spatial locations of the hauling locations (manure fields) 
sS2{2}=ExampleSources.ManureApp(:,1:2); 

% The spatial locations of the response (sampling sites) 
sR2=ExampleResponses.RABoBacsed(:,1:2);

% The response value 
zR2=ExampleResponses.RABoBacsed(:,3);

%% load information needed for the inform 
% The model type needs to be set to 3 for GORF
mt=3; 

% The type of proxy for flow we are using is Strahler Order flow, rather
% than upstream contributing river length
flowtype=2;

% We are not using the road network, instead using Euclidean distance proxy
useRoadNet=0;

% load the distance matrices and other needed information for the inform 
[D2,zS2,zRq2,ub2,lb2,modeltype2]=KewauneeObjFunLoad(sS2,zS2,sR2,zR2,flowtype,mt,useRoadNet);

%% get the parameter values 
% we inform independently of the meteorological factors 
STVs=0;

% we will use positive coefficients rather than R2 to inform
assumcase=1;

% we will set lower and upper bounds for the optimization based on the
% plots that were manually generated in the SI of the manuscript see 
% Figure S8(a) where we see more global optimums

% River Distance hyperparameter upper and lower bounds
Rbounds=[50 10000];

% Overland Flow hyper parameter upper and lower bounds
Obounds=[1000 18000];

% Ground Transport
Gbounds=[1000 18000];


% lower bounds
lb2 = [Gbounds(1), Obounds(1), Rbounds(1)];
ub2 = [Gbounds(2), Obounds(2), Rbounds(2)];

% We don't want to subset the data in any way
idxset=[]; 

% We want to use just pattern search to get the parameters, not bruteforce
bruteforce=0; 

% We want to use the penalty function to direct the parameter values away
% from things that have poor regression qualities 
usepenalty=1;

% INFORM!  
params2=SEDCparamOptim(D2,zS2,zRq2,STVs,assumcase,lb2,ub2,mt,idxset,bruteforce,usepenalty);

% We view the parameters that we found 

sprintf('gamma_G=%g km, alpha_O=%g km, alpha_R=%g m', params2(1)/1000, params2(2)/1000, params2(3))

%% obtain the HID source term values pre-z-scored

GORFafo2manureapp=KewauneeSEDCresults(D2,zS2,params2,zRq2,modeltype2);



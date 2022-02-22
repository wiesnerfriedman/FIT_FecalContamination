%% Load data
pathname=pwd; %the path where the file ARGdataPub was saved
q=sprintf('%s%s',pathname,'/FITBacteroides.mat'); 
load(q)

%% ORF example set sources!! and response data
% the source is high intensity developed land cover (HID)  
% and the response is the log10 relative abundance of human bacteroides 
% in sediment, but there are two coding options for HID. One is as gridded
% and equally weighted points. The other is as centroids of polygons with
% areas 

% CANDIDATE DATABASE 1- HID gridded points
% The spatial locations of the source
multisS1{1}=ExampleSources.HIDoptionA(:,1:2); 

% The scale of the source (here each equal to 1 because of being equally
% weighted gridded points)
multizS1{1}=ExampleSources.HIDoptionA(:,3); 

% CANDIDATE DATABASE 2- HID polygon centroids with polygon areas
% The spatial locations of the source
multisS1{2}=ExampleSources.HIDoptionB(:,1:2); 

% The scale of the source (here each equal to 1 because of being equally
% weighted gridded points)
multizS1{2}=ExampleSources.HIDoptionB(:,3); 

% The spatial locations of the response (sampling sites) 
sR1=ExampleResponses.RAHuBacSed(:,1:2);

% The response value 
zR1=ExampleResponses.RAHuBacSed(:,3);

%% Get set up for the Find 

% The model type needs to be set to 2 for ORF
mt=2; 

% The type of proxy for flow we are using is Strahler Order flow, rather
% than upstream contributing river length
flowtype=2;

% We are not using the road network, instead using Euclidean distance proxy
useRoadNet=0;

% River Distance hyperparameter upper and lower bounds
Rbounds=[249 3000]; 

% Overland Flow hyper parameter upper and lower bounds
Obounds=[100 3000];

% lower bounds
lb = [Obounds(1),Rbounds(1)];
ub = [Obounds(2),Rbounds(2)];

% We want to use just pattern search to get the parameters, not bruteforce
bruteforce=0; 

% We don't want to use the penalty function to direct the parameter values away
% from things that have poor regression qualities, because we do not know
% exactly what those would be
usepenalty=0;

% we will use positive coefficients rather than R2 to inform
assumcase=1;
%% Run the Find

RelInfo=ReliabilityScore(multisS1, multizS1, zR1, sR1, mt, flowtype, useRoadNet, assumcase,lb,ub, usepenalty, 6);


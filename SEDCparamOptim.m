function [params,ar,OFR]=SEDCparamOptim(D,zS,zRq,STVs,assumcase,lb,ub,mt,idxset,bruteforce,usepenalty)

% This function takes the distance matrix between sources and observations
% with associated observation flow & other values, as well as values
% associated with the scale, mass, or concentration at the source under the
% assumption of 0- no assumptions, 1- is a source, 2- somehting not valid
% anymore provided lower & upper bound values for the hyperparaemter values
% to be estimated. mt=1 is Euclidean, 2 is overland and river distance flow
%, 3 is ground transportation, overland and river distance flow. 
tic
% INPUT 
if nargin<10
    bruteforce=0; 
end 
if nargin<11
    usepenalty=0; 
end 

% D         distance matrix
% zS        source val
% zRq       observation vals
% assumcase source? or no assumption
% lb        hyperparameter lowerbound

A = [];
b = [];
Aeq = [];
beq = [];
options=optimoptions(@patternsearch,'MaxIterations',20,'Display','off');

% nsim=12;
% rng(1,'twister')% make results reproducible
% x0ops = lognrnd(0,1,nsim,1)*1000; % generate 8 lognormal options
% x0ops(x0ops>12000)=12000; % set maximum starting point of 12000
% x0ops(x0ops<50)=50; % set minimum starting point of 50
% rrx0{1}=x0ops; % for modeltype 1, use one dimensional random restart options
% 
% rng(3,'twister')% make results reproducible
% x0ops = lognrnd(0,1,nsim,2)*1000; % generate 8 lognormal options
% x0ops(x0ops>12000)=12000; % set maximum starting point of 12000
% x0ops(x0ops<50)=5; % set minimum starting point of 50
% rrx0{2}=x0ops; % for modeltype 2, use one dimensional random restart options
% 
% rng(5,'twister')% make results reproducible
% x0ops = lognrnd(0,1,nsim,3)*1000; % generate 8 lognormal options
% x0ops(x0ops>12000)=12000; % set maximum starting point of 12000
% x0ops(x0ops<50)=50; % set minimum starting point of 50
% rrx0{3}=x0ops; % for modeltype 2, use one dimensional random restart options

nsim=1;

zrq=zRq;

% get the values across the optimization=1, don't=0
%bruteforce=0;

switch mt
    case 1
        rng(1,'twister')% make results reproducible
        x0ops=randinterval(lb(1),ub(1),nsim);
        
        for qrr=1:size(x0ops,1) %%%%%%%%%%%% CAN THIS BE PARFOR?????
            x0=x0ops(qrr,1);
            [astar(qrr),fval(qrr)]=patternsearch(@(a) KewauneeObjectiveFunction1(D,zS,a,1,1,STVs,zrq,assumcase,idxset,1,usepenalty),x0,A,b,Aeq,beq,lb,ub,[],options);
        end
        
        %%% water %%%
        % find the best option between potential local
        % minima
        [~,idxmin]=min(fval);
        params=mean(astar(idxmin)); % make sure there is only one optimum if the resulting fval is the same between several parameter value
        
        % BRUTE FORCE METHOD
        if bruteforce==1
            ar=[1000:1000:10000]';
            for i=1:size(ar,1)
                [OFR(i),pval(i),R2(i)]=KewauneeObjectiveFunction1(D,zS,ar(i),1,1,STVs,zrq,assumcase,idxset,1,usepenalty);
            end
        end
        params=[params 1 1];
    case 2
        rng(1,'twister')% make results reproducible
        x0ops1=randinterval(lb(1),ub(1),nsim);
        x0ops2=randinterval(lb(2),ub(2),nsim);
        x0ops=[x0ops1 x0ops2]; % for modeltype 2, use one dimensional random restart options
        for qrr=1:size(x0ops,1) %%%%%%%%%%%% CAN THIS BE PARFOR?????
            x0=x0ops(qrr,:);
            [astar(qrr,:),fval(qrr)]=patternsearch(@(a) KewauneeObjectiveFunction1(D,zS,a(1),a(2),1,STVs,zrq,assumcase,idxset,3,usepenalty),x0,A,b,Aeq,beq,lb,ub,[],options);
        end
        %%% water %%%
        % find the best option between potential local
        % minima
        [~,idxmin]=min(fval);
        params(1)=mean(astar(idxmin,1)); % make sure there is only one optimum if the resulting fval is the same between several parameter values
        params(2)=mean(astar(idxmin,2));       
        
        params2(1)=5000; % make sure there is only one optimum if the resulting fval is the same between several parameter values
        params2(2)=3000;
        
        if bruteforce==1
            %ar=[0:0.1:1 2:249 250:5:12000 12005:20:24000]';
            ar=[50 1000:1000:10000]';
            for i=1:size(ar,1)
                [OFR(1,i),pval(1,i),R2(1,i)]=KewauneeObjectiveFunction1(D,zS,ar(i),params2(2),1,STVs,zrq,assumcase,idxset,3,usepenalty);
            end
            for i=1:size(ar,1)
                [OFR(2,i),pval(2,i),R2(2,i)]=KewauneeObjectiveFunction1(D,zS,params2(1),ar(i),1,STVs,zrq,assumcase,idxset,3,usepenalty);
            end
        end
        params=[params 1];
    case 3 
        rng(1,'twister')% make results reproducible
        x0ops1=randinterval(lb(1),ub(1),nsim);
        x0ops2=randinterval(lb(2),ub(2),nsim);
        x0ops3=randinterval(lb(3),ub(3),nsim);
        x0ops=[x0ops1 x0ops2 x0ops3]; % for modeltype 2, use one dimensional random restart options
        for qrr=1:size(x0ops,1) %%%%%%%%%%%% CAN THIS BE PARFOR?????
            x0=x0ops(qrr,:);
            [astar(qrr,:),fval(qrr)]=patternsearch(@(a) KewauneeObjectiveFunction1(D,zS,a(1),a(2),a(3),STVs,zrq,assumcase,idxset,5,usepenalty),x0,A,b,Aeq,beq,lb,ub,[],options);
        end
        %%% water %%%
        % find the best option between potential local
        % minima
        [~,idxmin]=min(fval);
        params(1)=mean(astar(idxmin,1)); % make sure there is only one optimum if the resulting fval is the same between several parameter values
        params(2)=mean(astar(idxmin,2));
        params(3)=mean(astar(idxmin,3));
        
        if bruteforce==1               
            ar=[0:0.1:1 2:249 250:5:12000 12005:20:24000]';
            for i=1:size(ar,1)
                [OFR(1,i),pval(1,i),R2(1,i)]=KewauneeObjectiveFunction1(D,zS,ar(i),params(2),params(3),STVs,zrq,assumcase,idxset,5,usepenalty);
            end
            for i=1:size(ar,1)
                [OFR(2,i),pval(2,i),R2(2,i)]=KewauneeObjectiveFunction1(D,zS,params(1),ar(i),params(3),STVs,zrq,assumcase,idxset,5,usepenalty);
            end
            for i=1:size(ar,1)
                [OFR(3,i),pval(3,i),R2(3,i)]=KewauneeObjectiveFunction1(D,zS,params(1),params(2),ar(i),STVs,zrq,assumcase,idxset,5,usepenalty);
            end
        end 
end 
toc 

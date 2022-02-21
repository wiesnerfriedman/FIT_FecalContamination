function [Dij,path1,d]=KewRoadDistanceMatrix(sS)

load('RoadNetworkMarch2021.mat','G','Len','NodeXY')

S1toNode=knnsearch(NodeXY,sS{1});
S2toNode=knnsearch(NodeXY,sS{2});

Dij=zeros(numel(S1toNode),numel(S2toNode));
for i=1:numel(S1toNode)
    n1=S1toNode(i); 
    for j=1:numel(S2toNode)
        n2=S2toNode(j); 
        [path1{i,j},d(i,j)] = shortestpath(G,n1,n2);
        Dij(i,j)=sum(Len(path1{i,j})); 
    end 
end 
function IDX=FindSitesonNetwork(RNpoints,sR,Qall)

DistMat=coord2dist(RNpoints,sR); 
for i=1:size(sR,1)
    tmp=find(DistMat(:,i)==0);
    if numel(tmp)>1
        [~,tmp2]=max(Qall(tmp));
        IDX(i)=tmp(tmp2); % if there are overlapping river reaches, picks the point with the greatest flow
    else
        IDX(i)=tmp;
    end
end
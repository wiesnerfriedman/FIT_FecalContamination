function result=GORFtinypartv2(Dkj,zS,gammaG,UseAk)


switch UseAk
    case 0
        mj=zS;
        
        A1=sedc3(Dkj,gammaG);
        A2=sum(A1,1);
        %A2(A2==0)=1;
        A=(A1./A2).*mj';
        result=sum(A,2);
    case 1
        
        mj=zS{1};
        Ak=zS{2};
        
        A1=sedc3(Dkj,gammaG);
        A2=sum(A1.*Ak,1);
        %A2(A2==0)=1;
        A=(A1./A2).*mj';
        result=sum(A,2).*Ak;
end

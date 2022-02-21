function D=sectionedcoord2dist(c1,c2)



n1=size(c1,1);
n2=size(c2,1);

if n1<=15000 && n2<=15000
    situation=1;
else
    situation=0;
end
    
    switch situation 
        case 1
            D=coord2dist(c1,c2);
            
        case 0
            thresh=10000;
            if n1>=n2
                incr=floor(n1/thresh);
                N=thresh*incr;
                D=[];
                remains=n1-N;
                for i=1:incr
                    p=(i*thresh-thresh)+1;
                    q=(i*thresh);
                    Dtemp=coord2dist(c1(p:q,:),c2);
                    D=[D;Dtemp];
                end
                Dtemp2=coord2dist(c1(end-remains+1:end,:),c2);
                D=[D;Dtemp2];
            end
            
            if n2>n1
                incr=floor(n2/thresh);
                N=thresh*incr;
                D=[];
                remains=n2-N;
                for i=1:incr
                    p=(i*thresh-thresh)+1;
                    q=(i*thresh);
                    Dtemp=coord2dist(c1,c2(p:q,:));
                    D=[D Dtemp];
                end
                Dtemp2=coord2dist(c1,c2(end-remains+1:end,:));
                D=[D Dtemp2];
            end
    end
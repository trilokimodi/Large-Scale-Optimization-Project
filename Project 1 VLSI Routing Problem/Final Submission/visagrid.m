function visagrid(dimx,dimy,nodelist,com,pi,shift)
    clf;
    hold on;
    node=1;
    for j=1:dimy
        for i=1:dimx
            if(pi(node)>1e-12)
                plot(i*100,100*j,'ro');
            else
                plot(i*100,100*j,'r.');
            end        
            node=node+1;    
        end
    end
    for i=1:dimx
        for j=1:dimy
            if(pi(node)>1e-12)
                plot(i*100+shift,100*j+shift,'bo');
            else
                plot(i*100+shift,100*j+shift,'b.');
            end    
            node=node+1;    
        end
    end
    axis([0 (dimx+1)*100 0 (dimy+1)*100]);
    
    last = 0;
    for dacom = 1 : size(com,1);
        first = last+1;
        slask = find(nodelist(last+1:length(nodelist)) == com(dacom,1));
        last = slask(1)+first-1;

        for i=first:last-1
            nodenr=nodelist(i);
            if nodenr<=dimx*dimy
                c1=100*[rem(nodenr-1,dimx)+1 floor((nodenr-1)/dimx)+1];
                c1=c1+[(dacom-1)*1 (dacom-1)*1];    
            else
                c1=100*[floor((nodenr-1-dimx*dimy)/dimy)+1 rem(nodenr-1-dimx*dimy,dimy)+1]; 
                c1=c1+[shift shift];    
                c1=c1+[(dacom-1)*1 (dacom-1)*1];    
            end    
            nodenr=nodelist(i+1);
            if nodenr<=dimx*dimy
                c2=100*[rem(nodenr-1,dimx)+1 floor((nodenr-1)/dimx)+1]; 
                c2=c2+[(dacom-1)*1 (dacom-1)*1];    
            else
                c2=100*[floor((nodenr-1-dimx*dimy)/dimy)+1 rem(nodenr-1-dimx*dimy,dimy)+1];
                c2=c2+[shift shift];   
                c2=c2+[(dacom-1)*1 (dacom-1)*1];    
            end    
            plot([c1(1) c2(1)],[c1(2) c2(2)]);            
        end
        
    end
end
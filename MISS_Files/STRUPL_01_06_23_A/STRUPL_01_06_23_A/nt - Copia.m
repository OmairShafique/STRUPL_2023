function [NT,ny]=nt(ELTYPE,ne,ragnetto)
if ragnetto==1
    ny=5;
    nrow=0;
    ncol=0;
    for i=1:ne
        if i==1
            nrow=nrow+1;
            NT(nrow,i)=1;
            
            nrow=nrow+1;
            NT(nrow,i)=-1;
            
        end
        if i==3
            nrow=nrow+1;
            NT(nrow,i)=1;
            
            nrow=nrow+1;
            NT(nrow,i)=-1;
            
        end
        if i==6
            nrow=nrow+1;
            NT(nrow,i)=1;
            
        end
    end
    
 

    
    
end
end

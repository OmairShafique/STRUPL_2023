

function [BS]=bsmx(ne,NDFPEREL,NSTRESSPEREL,ELTYPE)
%
    nrow=0;
    ncol=0;
    for i=1:ne
       nel=i;
       elt=ELTYPE(i);
       nste=NSTRESSPEREL(i);
       ndfe=NDFPEREL(i);
       [BE]=bemx(nel,nste,ndfe,elt);
%
        for k=1:nste
            nrow=nrow+1;
            for m=1:ndfe
                ncol=ncol+1;
                BS(nrow,ncol)=BE(k,m);
            end
        end
    end
    
        
            
        
         
       

end
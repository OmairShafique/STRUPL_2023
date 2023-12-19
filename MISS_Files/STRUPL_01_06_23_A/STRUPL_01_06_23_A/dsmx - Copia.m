function [DS]=dsmx(NSTRESSPEREL,ne,ELTYPE)
%
% computes matrix DSMX which collects, in diagonal form, the element
% elastic constitutive matrix DEMX times the element volume; it is not yet
% implemented the case of multile Gauss Points integration;
%
            nrow=0;
            ncol=0;
        for i=1:ne
            elt=ELTYPE(i);
            if elt>0
                return
            end
                nel=i;
                nste=NSTRESSPEREL(i);
                [DE]=demx(nel,nste,elt);
                    for k=1:nste
                        nrow=nrow+1;
                        for m=1:nste
                            ncol=ncol+1;
                            DS(nrow,ncol)=DE(k,m);
                            
                        end
                    end
         end
  end
        
        
     
        


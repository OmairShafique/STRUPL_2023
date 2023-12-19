

function [BS]=bsmx(ne,nn,ndf,COORD,NNODEPEREL,ELNODES,NSTRESSPEREL,ELTYPE)
%
% computes the diagnal element by element BS matrix of the structure
%
    nrow=0;
    ncol=0;
    for i=1:ne
       nel=i;
       elt=ELTYPE(i);
       nste=NSTRESSPEREL(i);
%
% Computes the BE matrix for element nel
%
           [BE,ndfe]=bemx(nel,nste,ndf,elt,COORD,ELNODES);
%
               for k=1:nste
                 nrow=nrow+1;
                 ncol=ndfe*(i-1);
                    for m=1:ndfe
                        ncol=ncol+1;
                        BS(nrow,ncol)=BE(k,m);
                    end
               end
      end
end
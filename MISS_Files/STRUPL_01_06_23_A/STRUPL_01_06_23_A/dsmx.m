function [DS,VOL]=dsmx(NSTRESSPEREL,ne,ELTYPE,ragnetto,...
            ELNODES,COORD,th)
%
% computes matrix DSMX which collects, in diagonal form, the element
% elastic constitutive matrix DEMX times the element volume; it is not yet
% implemented the case of multile Gauss Points integration;
%
            nrow=0;
            ncol=0;
        for i=1:ne
            elt=ELTYPE(i);
            nel=i;
            nste=NSTRESSPEREL(i);
            [DE]=demx(nel,nste,elt,ragnetto);
            if ragnetto==1
                VOL(i)=1;
                if nel==5
                    VOL(i)=2;
                end
                if nel==6
                    VOL(i)=2;
                end
            end
%            
%
           if elt==3
             [areatr]=areatriangle(nel,ELNODES,COORD);
             VOL(i)=areatr*th;
             DE=DE*VOL(i);
           end
%
               for k=1:nste
                 nrow=nrow+1;
                   for m=1:nste
                     ncol=nste*(i-1)+m;
                     DS(nrow,ncol)=DE(k,m);
                   end
               end
         end
  end
        
        
     
        


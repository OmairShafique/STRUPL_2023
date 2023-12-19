function [indexj,COLJ,modej]=inbasis(OUTB,ny,inrow,nno,...
    NOT,INB,flag3,indexjstandard,TA,LASTPAIR)
%
% display ('entro in inbasis')
%
% output: indexj identifies the pivot column of matrix TA
%
% inrow last mode which has exit from the Basis;
%
% 
% if inrow <=NY then incol=inrow+NY and indexj=k where k is such that
% OUTB(k)=incol
% if inrow >NY then incol= 2*NY-inrow and INDEXJ=k as before.
% 
% if the last activated mode coincides with the load factor alfa then
% indexj= 0 and the procedures terminates,
%
%
        for i=1:ny
            COLJ(i,1)=0;
        end
%
%
        if flag3==0
            for i=1:ny+1
                if OUTB(i)==LASTPAIR(1,2)
                    indexj=i+1
                    modej=OUTB(i)
                end % if OUTB
            end% for i
         end %if flag3
%
 if flag3==1
     jmode=OUTB(ny+1)-1;
     modej=jmode;
     for k=1:ny+1
        if OUTB(k)==jmode
            indexj=k+1;
        end
     end
 end
%
% Computes the column of TA corresponding to the variable entering into the
% Basis
%
   for k=1:ny
       COLJ(k,1)=TA(k,indexj);
    end
%
% display ('esco da inbasis')
%
% indexj=indexj
end

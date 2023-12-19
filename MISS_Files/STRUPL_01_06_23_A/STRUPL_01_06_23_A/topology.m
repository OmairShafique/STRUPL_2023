function [TS]=topology(ndf,nn,ELNODES,ne,maxnnperel,...
        NNODEPEREL,ELTYPE)
%
% 
% Computes matrix TS topology of the structure which transforms elemnt 
% displacements into nodal dispalcements;
% NNODEPEREL lists the number of nodes per elemnt;
%
    rownumber=0;
        for i=1:ne
        numnode=NNODEPEREL(i);
            for j=1:numnode
                nodenumber=ELNODES(i,j);
                 nodenm1=nodenumber-1;
                    tndfm1=0;
                    if nodenm1>0
                        for m=1:nodenm1
                          tndfm1=tndfm1+ndf;                       
                        end
                    end      
                            
                         for k=1:ndf
                            rownumber=rownumber+1;
                            columnnumber=tndfm1+k;
                            TS(rownumber,columnnumber)=1;
                         end
               end
        end 
 
    end  
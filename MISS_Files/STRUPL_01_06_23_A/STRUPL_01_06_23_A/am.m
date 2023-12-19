function [AS,maxasii]=am(NT,DS,BS,TS,TF,KSFMM1,TFM,HS,RCR,TCR,ragnetto, ...
    icrk,KE,ny)
%
%   computes matrix AS
%
    if ragnetto==1
        AS=HS+NT*(DS-DS*BS*TS*TF*TFM*KSFMM1*TFM'*TF'*TS'*BS'*DS)*NT';
    end
    if icrk==1
%
     AS=HS+NT*(RCR*TCR*KE*TCR'*RCR'-RCR*TCR*KE*TS*TF*TFM*KSFMM1*TFM'*TF'*TS'...
    *KE'*TCR'*RCR')*NT';
%
    end %if icrk
%    
% definition of the max diagonal term
%
        maxasii=0
        for i =1:ny
            asii=AS(i,i);
            if asii>maxasii
                maxasii=asii;
            end % if
        end %for i
   
end
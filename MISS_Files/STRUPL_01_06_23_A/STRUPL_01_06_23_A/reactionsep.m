
function [RSEP]=reactionsep(UMEP,UC,TC,TF,KS,TFM,TS,BS,DS,NT,LAMBDA,...
    icrk,KE,TCR,RCR,FS)
%
% computes nodal reactionsdue to imposed displacements at the constrained
% dof and nodal displacement at the free dof
%
%
   
   
   if icrk==0
    RSEP=TC'*KS*TF*TFM*UMEP+TC'*KS*TC*UC-TC'*TS'*BS'*DS*NT'*LAMBDA-TC'*FS(:,1);
   end
   if icrk==1
    RSEP=TC'*KS*UMEP+TC'*KS*TC*UC-TC'*TS'*KE'*TCR'*RCR'*NT'*LAMBDA-TC'*FS(:,1);
   end 
   
   
            
%
%
%
end

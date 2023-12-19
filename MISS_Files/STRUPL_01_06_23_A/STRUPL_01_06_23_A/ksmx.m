
function [KSFM,KSF,KS,KE]=ksmx(TS,TF,BS,DS,TFM)
%
% computes structure stiffnes matrix KE relative to all elements
%
    KE=BS'*DS*BS;
%
% computes structure stiffness matrix relative to all dof;
%
    KS=TS'*KE*TS; 
%
% computes structural stiffness matrix relative to the free dof;
%
    KSF=TF'*KS*TF;
    
%
% computes structural stiffness matrix relative to master free dof
%
    KSFM=TFM'*KSF*TFM;
%
%        
        end
function [BI]=bi(ny,NT,DS,BS,TS,TF,TFM,KSM1,FS,ragnetto,icrk,KE,...
TCR,RCR)
%
% Computes vector BI=N'Qel
%
    if ragnetto==1
        BI=NT*DS*BS*TS*TF*TFM*KSM1*FS;
    end
    if icrk==1
      BI=NT*RCR*TCR*KE*TS*TF*TFM*KSM1*FS;
    end  
end

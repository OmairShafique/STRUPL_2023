function [BI]=bi(ny,NT,DS,BS,TS,VS,KSM1,FS)
%
% Computes vector BI=N'Qel
%
    BI=NT*DS*BS*TS*VS*KSM1*FS;
end

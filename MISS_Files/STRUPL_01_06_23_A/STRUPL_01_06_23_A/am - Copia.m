function [AS]=am(NT,DS,BS,TS,VS,KSM1)
%
%   computes matrix AS
%
    AS=NT*(DS-DS*BS*TS*VS*KSM1*VS'*TS'*BS'*DS)*NT';
    
    
end
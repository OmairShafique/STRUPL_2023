        function [KS]=ksmx(TS,VS,BS,DS)
%
% compute structure stiffnes matrix KS
%
        KS=VS'*TS'*BS'*DS*BS*TS*VS;

%
%
        end
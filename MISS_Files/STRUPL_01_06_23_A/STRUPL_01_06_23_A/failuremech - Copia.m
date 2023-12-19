function [YDOT]=failuremech(ny,TA,indexj);
%
%test for collapse mechanism
% 
%
% display ('In failuremech')

%
%

% display ('collapse mechanism')
                     nym1=ny-1
                     
                     YDOT(ny,1)=1;
                     YDOT(1:nym1,1)=TA(1:nym1,indexj);
% display ('out failuremech')
                     

             end
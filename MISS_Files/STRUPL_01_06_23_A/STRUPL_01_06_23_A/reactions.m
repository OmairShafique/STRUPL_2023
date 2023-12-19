function [RS]=reactions(US,UC,TC,TF,KS,TFM,FS)
%
% computes nodal reactions due to imposed displacements at the constrained
% dof and nodal displacement at the free dof and due to loads applied
% directely to the supports;
%
% V1=TFM*US
% V2=TF*V1
% V3=KS*V2
% V4=TC'*V3
% V5=TC*UC
% V6=KS*V5
% V7=TC'*V6
                RS=(TC')*KS*TF*TFM*US+(TC')*KS*TC*UC+(TC')*FS(:,1);
%
%
%
end

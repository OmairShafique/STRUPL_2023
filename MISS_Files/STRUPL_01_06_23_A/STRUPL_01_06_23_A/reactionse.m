function [RSE]=reactionse(UME,UC,TC,TF,KS,TFM,FS)
%
% computes nodal reactionsdue to imposed displacements at the constrained
% dof and nodal displacement at the free dof
%
%
                RCC=TC'*FS(:,1)
                RSE=(TC')*KS*UME-RCC
%
%
%
end

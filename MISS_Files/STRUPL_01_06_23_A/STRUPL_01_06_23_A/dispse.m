
function [UME]= dispse (KSFMM1,FSFM,TF,KS,TC,UC,TFM)
%
%
%   computes the free master elastic displacements UM due to the free master nodal 
%   loads FSFMand given prescribed nodal prescribed displacements UC 
%
%
%    
    UME=TF*TFM*KSFMM1*(FSFM-TFM'*TF'*KS*TC*UC)+TC*UC;
%
%
end
    
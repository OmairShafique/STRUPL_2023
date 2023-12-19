function [US]= disps (KSFMM1,FSFM,TF,KS,TC,UC,TFM)
%
%
%   computes the free displacements US due to the free nodal loads FSF,
%   given the nodal loads FS relative to all structural nodes
%
%
%    
    US=KSFMM1*(FSFM-TFM'*TF'*KS*TC*UC);
%
%
end
    
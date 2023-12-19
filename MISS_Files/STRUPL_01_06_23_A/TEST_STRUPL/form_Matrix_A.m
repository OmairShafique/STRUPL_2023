function [A_Matrix]=form_Matrix_A(Normality_Matrix,dee,bee,KK)
%
% computes matrix AS

%
 A_Matrix=Hard_Soft_Matrix+Normality_Matrix*(dee-dee*bee*inv(KK)*bee'*dee)*Normality_Matrix';


%  A_Matrix=Hard_Soft_Matrix+Normality_Matrix*(Rotation_Matrix*TCR*kk*TCR'*Rotation_Matrix'-Rotation_Matrix*TCR*KE*TS*TF*TFM*kk*TFM'*TF'*TS'...
  %  *kk'*TCR'*Rotation_Matrix')*Normality_Matrix';

% AS=HS+NT*(RCR*TCR*KE*TCR'*RCR'-RCR*TCR*KE*TS*TF*TFM*KSFMM1*TFM'*TF'*TS'*KE'*TCR'*RCR')*NT';

% definition of the max diagonal term
%
        maxasii=0
        for i =1:ny
            asii=AS(i,i);
            if asii>maxasii
                maxasii=asii;
            end % if
        end %for i
   
end
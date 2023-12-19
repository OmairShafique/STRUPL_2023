function [ROT] = Rotation_matrix( beta, v )
%Calcola la matrice di rotazione del sistema di riferimento locale n,t rispetto
%al sistema di riferimento globale x1,x2 (verso di rotazione di x1 su x2 sempre antiorario).
%-------------------------------------------------------------
% INPUT:    
%           beta   :  angolo tra asse x1 e n
%           v      :  verso di rotazione di n su t 
%                     1 --> orario
%                     2 --> anti-orario
% 
% OUTPUT:   ROT    :  Matrice di rotazione che permette di passare dal
%                     sistema di riferimento globale a locale

%--------------------------------------------------------------

if v == 1
    ROT = [cos(beta) sin(beta);
            sin(beta) -cos(beta)];
elseif v == 2
    ROT = [cos(beta) sin(beta);
            -sin(beta) cos(beta)];
end
    
end


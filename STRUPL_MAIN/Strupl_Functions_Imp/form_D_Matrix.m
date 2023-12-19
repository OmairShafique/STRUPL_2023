function [D_Matrix]=form_D_Matrix (Elastic_Modulus,Poisson_Ratio,thickness_of_Plate,Element_type,ngpb,ngps)
%
% This function forms the elasticity matrix in a plate element
% concerning also shear and bending actions
%
if Element_type==3 && ngpb==0
    c=Elastic_Modulus/(1.-Poisson_Ratio*Poisson_Ratio);
    %
    D_Matrix=c*[1 Poisson_Ratio 0. ;...
        Poisson_Ratio 1 0. ;...
        0. 0. .5*(1.-Poisson_Ratio)];

elseif Element_type~=3 && ngpb~=0
    DR= Elastic_Modulus*(thickness_of_Plate^3)/(12*(1.-Poisson_Ratio*Poisson_Ratio));
    %
    D_Matrix=DR*[1 Poisson_Ratio 0. ;...
        Poisson_Ratio 1 0. ;...
        0. 0. (1.-Poisson_Ratio)/2] ;
elseif Element_type~=3 && ngps~=0
    G= Elatic_Modulus/(2*(1.+Poisson_Ratio));
    %
    D_Matrix=G* [thickness_of_Plate     0 ;...
        0     thickness_of_Plate];
end
% end function form_D_Matrix
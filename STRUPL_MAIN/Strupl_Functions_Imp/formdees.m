function[dees] = formdees(Elatic_Modulus,Poisson_Ratio,thickness_of_Plate)
%
% This function forms the elasticity matrix for the shear
% action in a thick plate element
%
G= Elatic_Modulus/(2*(1.+Poisson_Ratio));
%
dees=G* [thickness_of_Plate     0 ;...
          0     thickness_of_Plate];
end
function[deeb] = formdeeb(Elatic_Modulus,Poisson_Ratio,thickness_of_Plate)
%
% This function forms the elasticity matrix for a bending
% action in a plate element
%
DR= Elatic_Modulus*(thickness_of_Plate^3)/(12*(1.-Poisson_Ratio*Poisson_Ratio));
%
deeb=DR*[1 Poisson_Ratio 0. ;...
Poisson_Ratio 1 0. ;...
0. 0. (1.-Poisson_Ratio)/2] ;
%
% end function fromdeeb
end
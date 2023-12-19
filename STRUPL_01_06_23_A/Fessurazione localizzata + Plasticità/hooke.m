function [D]=hooke(at,E,nu)
% HOOKE --> Calculate the material matrix for a linear
% elastic and isotropic material.
%
% INPUT:    ptype=1 : Plane stress
%               2   : Plane strain
%               3   : Axisymmetry
%               4   : Three dimensional
%
%           E  : Young's modulus
%           v  : Poissons const.
%
% OUTPUT:   D  : Material matrix
%-------------------------------------------------------------

 if at==1
        Dm=E/(1-nu^2)*[1  nu   0;
                      nu  1   0;
                      0  0 (1-nu)/2];
 elseif at==2
        Dm=E/(1+nu)/(1-2*nu)*[1-nu  nu    nu        0;
                             nu  1-nu   nu        0;
                             nu   nu   1-nu       0;
                             0   0    0   (1-2*nu)/2];;
 elseif at==3
        Dm=E/(1+nu)/(1-2*nu)*[1-nu  nu    nu        0;
                             nu  1-nu   nu        0;
                             nu   nu   1-nu       0;
                             0   0    0   (1-2*nu)/2];;
 elseif at==4
        Dm=E/(1+nu)/(1-2*nu)*[1-nu  nu    nu    0    0    0;
                             nu  1-nu   nu    0    0    0;
                             nu   nu   1-nu   0    0    0;
                             0   0    0 (1-2*nu)/2    0    0;
                             0   0    0    0 (1-2*nu)/2    0;
                             0   0    0    0    0  (1-2*nu)/2];
 else
   error('Error ! Check first argument, ptype=1,2,3 or 4 allowed')
   return
 end
 D=Dm;
%--------------------------end--------------------------------


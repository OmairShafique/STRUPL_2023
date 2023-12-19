function [Ke,fe]=stif2d(et,ex,ey,t,b,D)
% STIF2D --> The routine computes the stiffness matrix for a 3 node  or 
%  4 nodes isoparametric element in plane stress 
%-------------------------------------------------------------
% INPUT:    at : Analysis type 1=plane stress  2=plane strain
%           et : Element type: 1=4 node quad 2=3 node triangle
%           ex : Nodes coordinates in x direction
%           ey : Nodes coordinates  in y direction
%           t  : Element thickness
%           b  : Body force vector in x and y direction 
%           D  : Constitutive matrix
%
% OUTPUT:   Ke : Element stiffness matrix
%           fe : Equivalent nodal forces
%-------------------------------------------------------------

tol = 10^(-8) ;

if et == 1 % 4 nodes quad
    nn = 4 ;
    % Gauss points coordinates and weights(ngp=4)-------------
    ngp = 4 ;
    r = 0.577350269189626 ; w1 = 1 ;
    gp(:,1)=[-r; r; -r; r];  gp(:,2)=[-r; -r; r; r];
    w(:,1)=[ w1; w1; w1; w1];   w(:,2)=[ w1; w1; w1; w1];
    wp=w(:,1).*w(:,2);
    xsi=gp(:,1);  eta=gp(:,2);    
    % Shape function -----------------------------------------
    N(:,1)=(1-xsi).*(1-eta)/4;  N(:,3)=(1+xsi).*(1+eta)/4;  
    N(:,2)=(1+xsi).*(1-eta)/4;  N(:,4)=(1-xsi).*(1+eta)/4;
    % Shape function derivative ------------------------------
    %respect to xsi
    dNr(1:2:2*ngp,1)=-(1-eta)/4; dNr(1:2:2*ngp,2)= (1-eta)/4;  
    dNr(1:2:2*ngp,3)= (1+eta)/4; dNr(1:2:2*ngp,4)= -(1+eta)/4;
    %respect to eta
    dNr(2:2:2*ngp,1)=-(1-xsi)/4; dNr(2:2:2*ngp,2)=-(1+xsi)/4;  
    dNr(2:2:2*ngp,3)= (1+xsi)/4; dNr(2:2:2*ngp,4)= (1-xsi)/4;      
elseif et == 2 % 3 nodes triangle
    nn = 3 ;
    % Gauss points coordinates and weights(ngp=1)-------------
    ngp = 1 ; 
    r = 0.333333333333333 ; w = 0.5 ; wp = w ;
    gp = [r r]; 
    xsi=gp(:,1);  eta=gp(:,2);  
    % Shape function -----------------------------------------
    N(:,1)= 1 - xsi - eta; 
    N(:,2)= xsi;
    N(:,3)= eta;  
    % Shape function derivative ------------------------------
    %respect to xsi       respect to eta
    dNr(1:2:2*ngp,1)= -1; dNr(2:2:2*ngp,1)= -1; 
    dNr(1:2:2*ngp,2)= 1;  dNr(2:2:2*ngp,2)= 0; 
    dNr(1:2:2*ngp,3)= 0;  dNr(2:2:2*ngp,3)= 1;   
end

JT = dNr*[ex;ey]' ; % Complete Jacobian matrix
Ke=zeros(2*nn);
fe=zeros(2*nn,1);

% Computation of B compatibility matrix  ---------------------
  for i=1:ngp
      indx=[ 2*i-1; 2*i ];
      detJ=det(JT(indx,:)) ;
      if detJ<tol
        disp('Jacobideterminant equal or less than zero!')
      end
      JTinv=inv(JT(indx,:));
      dNx=JTinv*dNr(indx,:);
      
      B(1,1:2:2*nn-1)=dNx(1,:);
      B(2,2:2:2*nn)  =dNx(2,:);
      B(3,1:2:2*nn-1)=dNx(2,:);
      B(3,2:2:2*nn)  =dNx(1,:);

      N2(1,1:2:2*nn-1)=N(i,:) ;
      N2(2,2:2:2*nn)  =N(i,:) ;
      
      Ke=Ke+B'*D*B*detJ*wp(i)*t ;
      fe=fe+N2'*b*detJ*wp(i)*t ;
  end
%--------------------------end--------------------------------


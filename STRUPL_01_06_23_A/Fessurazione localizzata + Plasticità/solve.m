  function [u,R]=solve(K,F,bc)
% SOLVE --> Solve static FE-equations considering boundary conditions.
%-------------------------------------------------------------
% INPUT:   K  : global stiffness matrix, dim(K)= nd x nd, nd : number of dof's
%          F  : global load vector, dim(f)= nd x 1
%          bc : boundary condition matrix
%               dim(bc)= nbc x 2, nbc = number of b.c.'s
%               bc = (gdl pd), pd = prescribed displacemnt
%
% OUTPUT:  u : displacement solution including boundary values
%          R : reaction force vector
%              dim(u)=dim(R)= nd x 1 
%-------------------------------------------------------------
  if nargin==2 ; 
     u=K\F;  
  elseif nargin==3;
     [nd,nd]=size(K);
     fdof=[1:nd]';  %free dof
     
     u=zeros(size(fdof));
     
     pdof=bc(:,1);  %prescribed dof
     pd=bc(:,2);    %prescribed displacement
     fdof(pdof)=[]; 

     s=K(fdof,fdof)\(F(fdof)-K(fdof,pdof)*pd);
    
     u(pdof)=pd;
     u(fdof)=s;
  end  
     R=K*u-F;
%--------------------------end--------------------------------
